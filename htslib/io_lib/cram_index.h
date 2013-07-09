#ifndef _CRAM_INDEX_H_
#define _CRAM_INDEX_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Loads a CRAM .crai index into memory.
 * Returns 0 for success
 *        -1 for failure
 */
int cram_index_load(cram_fd *fd, char *fn);

void cram_index_free(cram_fd *fd);

/*
 * Searches the index for the first slice overlapping a reference ID
 * and position.
 *
 * Returns the cram_index pointer on sucess
 *         NULL on failure
 */
cram_index *cram_index_query(cram_fd *fd, int refid, int pos, cram_index *frm);

/*
 * Skips to a container overlapping the start coordinate listed in
 * cram_range.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek_to_refpos(cram_fd *fd, cram_range *r);

/*
 * Seek within a cram file.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int cram_seek(cram_fd *fd, off_t offset, int whence);

#ifdef __cplusplus
}
#endif

#endif
