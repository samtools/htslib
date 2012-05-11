#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include "endian.h"
#include "bgzf.h"
#include "sam.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

int sam_verbose = 3;

/******************
 * Basic file I/O *
 ******************/

samFile *sam_open(const char *fn, const char *mode, const char *fn_ref)
{ // nearly identical to vcf_open()
	const char *p;
	samFile *fp;
	fp = (samFile*)calloc(1, sizeof(samFile));
	for (p = mode; *p; ++p) {
		if (*p == 'w') fp->is_write = 1;
		else if (*p == 'b') fp->is_bin = 1;
	}
	if (fp->is_bin) {
		if (fp->is_write) fp->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_dopen(fileno(stdout), mode);
		else fp->fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
	} else {
		if (fp->is_write) {
			fp->fp = strcmp(fn, "-")? fopen(fn, "rb") : stdout;
		} else {
			gzFile gzfp;
			gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
			if (gzfp) fp->fp = ks_init(gzfp);
			if (fn_ref) fp->fn_ref = strdup(fn_ref);
		}
	}
	if (fp->fp == 0) {
		if (sam_verbose >= 2)
			fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);
		free(fp->fn_ref); free(fp);
		return 0;
	}
	return fp;
}

void sam_close(samFile *fp)
{
	if (!fp->is_bin) {
		free(fp->line.s);
		if (!fp->is_write) {
			gzFile gzfp = ((kstream_t*)fp->fp)->f;
			ks_destroy((kstream_t*)fp->fp);
			gzclose(gzfp);
			free(fp->fn_ref);
		} else fclose((FILE*)fp->fp);
	} else bgzf_close((BGZF*)fp->fp);
	free(fp);
}

/**************
 * Header I/O *
 **************/

sam_hdr_t *sam_hdr_init()
{
	sam_hdr_t *h;
	h = (sam_hdr_t*)calloc(1, sizeof(sam_hdr_t));
	h->is_be = ed_is_big();
	return h;
}

void sam_hdr_destroy(sam_hdr_t *h)
{
}

sam_hdr_t *sam_hdr_read(samFile *fp)
{
	sam_hdr_t *h;
	if (fp->is_bin) {
		char buf[4];
		int magic_len, has_EOF;
		int32_t i = 1, name_len;
		// check EOF
		has_EOF = bgzf_check_EOF((BGZF*)fp->fp);
		if (has_EOF < 0) {
			if (errno != ESPIPE && sam_verbose >= 2) perror("[W::sam_hdr_read] bgzf_check_EOF");
		} else if (has_EOF == 0 && sam_verbose >= 2)
			fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
		// read "BAM1"
		magic_len = bgzf_read((BGZF*)fp->fp, buf, 4);
		if (magic_len != 4 || strncmp(buf, "BAM\1", 4)) {
			if (sam_verbose >= 1) fprintf(stderr, "[E::%s] invalid BAM binary header\n", __func__);
			return 0;
		}
		h = sam_hdr_init();
		// read plain text and the number of reference sequences
		bgzf_read((BGZF*)fp->fp, &h->l_text, 4);
		if (h->is_be) ed_swap_4p(&h->l_text);
		h->text = (char*)calloc(h->l_text + 1, 1);
		bgzf_read((BGZF*)fp->fp, h->text, h->l_text);
		bgzf_read((BGZF*)fp->fp, &h->n_targets, 4);
		if (h->is_be) ed_swap_4p(&h->n_targets);
		// read reference sequence names and lengths
		h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
		h->target_len = (uint32_t*)calloc(h->n_targets, 4);
		for (i = 0; i != h->n_targets; ++i) {
			bgzf_read((BGZF*)fp->fp, &name_len, 4);
			if (h->is_be) ed_swap_4p(&name_len);
			h->target_name[i] = (char*)calloc(name_len, 1);
			bgzf_read((BGZF*)fp->fp, h->target_name[i], name_len);
			bgzf_read((BGZF*)fp->fp, &h->target_len[i], 4);
			if (h->is_be) ed_swap_4p(&h->target_len[i]);
		}
		return h;
	}
	return 0;
}
