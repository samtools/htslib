/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 *
 * Support for CRAM index format: foo.cram.crai
 */

/*
 * The index is a gzipped tab-delimited text file with one line per slice.
 * The columns are:
 * 1: reference number (0 to N-1, as per BAM ref_id)
 * 2: reference position of 1st read in slice (1..?)
 * 3: number of reads in slice
 * 4: offset of container start (relative to end of SAM header, so 1st
 *    container is offset 0).
 * 5: slice number within container (ie which landmark).
 *
 * In memory, we hold this in a nested containment list. Each list element is
 * a cram_index struct. Each element in turn can contain its own list of
 * cram_index structs.
 *
 * Any start..end range which is entirely contained within another (and
 * earlier as it is sorted) range will be held within it. This ensures that
 * the outer list will never have containments and we can safely do a
 * binary search to find the first range which overlaps any given coordinate.
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"
#include "io_lib/zfio.h"

#if 0
static void dump_index_(cram_index *e, int level) {
    int i, n;
    n = printf("%*s%d / %d .. %d, ", level*4, "", e->refid, e->start, e->end);
    printf("%*soffset %"PRId64"\n", MAX(0,50-n), "", e->offset);
    for (i = 0; i < e->nslice; i++) {
	dump_index_(&e->e[i], level+1);
    }
}

static void dump_index(cram_fd *fd) {
    int i;
    for (i = 0; i < fd->index_sz; i++) {
	dump_index_(&fd->index[i], 0);
    }
}
#endif

/*
 * Loads a CRAM .crai index into memory.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, char *fn) {
    char line[1024], fn2[PATH_MAX];
    zfp *fp;
    cram_index *idx;
    cram_index **idx_stack = NULL, *ep, e;
    int idx_stack_alloc = 0, idx_stack_ptr = 0;

    /* Check if already loaded */
    if (fd->index)
	return 0;

    fd->index = calloc((fd->index_sz = 1), sizeof(*fd->index));
    if (!fd->index)
	return -1;

    idx = &fd->index[0];
    idx->refid = -1;
    idx->start = INT_MIN;
    idx->end   = INT_MAX;

    idx_stack = calloc(++idx_stack_alloc, sizeof(*idx_stack));
    idx_stack[idx_stack_ptr] = idx;

    sprintf(fn2, "%s.crai", fn);
    if (!(fp = zfopen(fn2, "r"))) {
	perror(fn2);
	return -1; 
    }

    while (zfgets(line, 1024, fp)) {
	/* 1.1 layout */
	sscanf(line, "%d\t%d\t%d\t%"PRId64"\t%d\t%d",
	       &e.refid,
	       &e.start,
	       &e.end,
	       &e.offset,
	       &e.slice,
	       &e.len);
	e.end += e.start-1;
	//printf("%d/%d..%d\n", e.refid, e.start, e.end);

	if (e.refid != idx->refid) {
	    if (fd->index_sz < e.refid+2) {
		fd->index_sz = e.refid+2;
		fd->index = realloc(fd->index,
				    fd->index_sz * sizeof(*fd->index));
	    }
	    idx = &fd->index[e.refid+1];
	    idx->refid = e.refid;
	    idx->start = INT_MIN;
	    idx->end   = INT_MAX;
	    idx->nslice = idx->nalloc = 0;
	    idx->e = NULL;
	    idx_stack[(idx_stack_ptr = 0)] = idx;
	}

	while (!(e.start >= idx->start && e.end <= idx->end)) {
	    idx = idx_stack[--idx_stack_ptr];
	}

	// Now contains, so append
	if (idx->nslice+1 >= idx->nalloc) {
	    idx->nalloc = idx->nalloc ? idx->nalloc*2 : 16;
	    idx->e = realloc(idx->e, idx->nalloc * sizeof(*idx->e));
	}

	e.nalloc = e.nslice = 0; e.e = NULL;
	*(ep = &idx->e[idx->nslice++]) = e;
	idx = ep;

	if (++idx_stack_ptr >= idx_stack_alloc) {
	    idx_stack_alloc *= 2;
	    idx_stack = realloc(idx_stack, idx_stack_alloc*sizeof(*idx_stack));
	}
	idx_stack[idx_stack_ptr] = idx;
    }
    zfclose(fp);
    free(idx_stack);

    // dump_index(fd);

    return 0;
}

static void cram_index_free_recurse(cram_index *e) {
    if (e->e) {
	int i;
	for (i = 0; i < e->nslice; i++) {
	    cram_index_free_recurse(&e->e[i]);
	}
	free(e->e);
    }
}

void cram_index_free(cram_fd *fd) {
    int i;

    if (!fd->index)
	return;
    
    for (i = 0; i < fd->index_sz; i++) {
	cram_index_free_recurse(&fd->index[i]);
    }
    free(fd->index);

    fd->index = NULL;
}

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position, or one immediately preceeding it if none is found in
 * the index to overlap this position. (Our index may have missing
 * entries, but we require at least one per reference.)
 *
 * If the index finds multiple slices overlapping this position we
 * return the first one only. Subsequent calls should specifying
 * "from" as the last slice we checked to find the next one. Otherwise
 * set "from" to be NULL to find the first one.
 *
 * Returns the cram_index pointer on sucess
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, int pos, 
			     cram_index *from) {
    int i, j, k;
    cram_index *e;

    if (refid+1 < 0 || refid+1 >= fd->index_sz)
	return NULL;

    i = 0, j = fd->index[refid+1].nslice-1;

    if (!from)
	from = &fd->index[refid+1];

    for (k = j/2; k != i; k = (j-i)/2 + i) {
	if (from->e[k].refid > refid) {
	    j = k;
	    continue;
	}

	if (from->e[k].refid < refid) {
	    i = k;
	    continue;
	}

	if (from->e[k].start >= pos) {
	    j = k;
	    continue;
	}

	if (from->e[k].start < pos) {
	    i = k;
	    continue;
	}
    }

    /* The above found *a* bin overlapping, but not necessarily the first */
    while (i > 0 && from->e[i-1].end >= pos)
	i--;

    /* Special case for matching a start pos */
    if (i+1 < from->nslice &&
	from->e[i+1].start == pos &&
	from->e[i+1].refid == refid)
	i++;

    e = &from->e[i];

    return e;
}

/*
 * Seek within a cram file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence) {
    char buf[65536];

    if (fseeko(fd->fp, offset, whence) == 0)
	return 0;

    if (!(whence == SEEK_CUR && offset >= 0))
	return -1;

    /* Couldn't fseek, but we're in SEEK_CUR mode so read instead */
    while (offset > 0) {
	int len = MIN(65536, offset);
	if (len != fread(buf, 1, len, fd->fp))
	    return -1;
	offset -= len;
    }

    return 0;
}


/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * In theory we call cram_index_query multiple times, once per slice
 * overlapping the range. However slices may be absent from the index
 * which makes this problematic. Instead we find the left-most slice
 * and then read from then on, skipping decoding of slices and/or
 * whole containers when they don't overlap the specified cram_range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r) {
    cram_index *e;

    // Ideally use an index, so see if we have one.
    if ((e = cram_index_query(fd, r->refid, r->start, NULL))) {
	if (0 != cram_seek(fd, e->offset, SEEK_SET))
	    if (0 != cram_seek(fd, e->offset - fd->first_container, SEEK_CUR))
		return -1;
    } else {
	fprintf(stderr, "Unknown reference ID. Missing from index?\n");
	return -1;
    }

    if (fd->ctr) {
	cram_free_container(fd->ctr);
	fd->ctr = NULL;
    }

    return 0;
}
