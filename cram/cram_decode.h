/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the decoding portion of CRAM I/O. Also see
 * cram_codecs.[ch] for the actual encoding functions themselves.
 */

#ifndef _CRAM_READ_H_
#define _CRAM_READ_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * CRAM sequence iterators.
 */

/*! Read the next cram record and return it as a cram_record.
 *
 * Note that to decode cram_record the caller will need to look up some data
 * in the current slice, pointed to by fd->ctr->slice. This is valid until
 * the next call to cram_get_seq (which may invalidate it).
 *
 * @return
 * Returns record pointer on success (do not free);
 *        NULL on failure
 */
cram_record *cram_get_seq(cram_fd *fd);

/*! Read the next cram record and convert it to a bam_seq_t struct.
 *
 * @return
 * Returns 0 on success;
 *        -1 on EOF or failure (check fd->err)
 */
int cram_get_bam_seq(cram_fd *fd, bam_seq_t **bam);


/* ----------------------------------------------------------------------
 * Internal functions
 */

/*! INTERNAL:
 * Decodes a CRAM block compression header.
 *
 * @return
 * Returns header ptr on success;
 *         NULL on failure
 */
cram_block_compression_hdr *cram_decode_compression_header(cram_fd *fd,
							   cram_block *b);

/*! INTERNAL:
 * Decodes a CRAM (un)mapped slice header block.
 *
 * @return
 * Returns slice header ptr on success;
 *         NULL on failure
 */
cram_block_slice_hdr *cram_decode_slice_header(cram_fd *fd, cram_block *b);


/*! INTERNAL:
 * Decode an entire slice from container blocks. Fills out s->crecs[] array.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
		      SAM_hdr *hdr);


#ifdef __cplusplus
}
#endif

#endif
