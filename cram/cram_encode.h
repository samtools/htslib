/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the encoding portion of CRAM I/O. Also see
 * cram_codecs.[ch] for the actual encoding functions themselves.
 */

#ifndef _CRAM_WRITE_H_
#define _CRAM_WRITE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * CRAM sequence iterators.
 */

/*! Write iterator: put BAM format sequences into a CRAM file.
 *
 * We buffer up a containers worth of data at a time.
 *
 * FIXME: break this into smaller pieces.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_put_bam_seq(cram_fd *fd, bam_seq_t *b);


/* ----------------------------------------------------------------------
 * Internal functions
 */

/*! INTERNAL:
 * Encodes a compression header block into a generic cram_block structure.
 *
 * @return
 * Returns cram_block ptr on success;
 *         NULL on failure
 */
cram_block *cram_encode_compression_header(cram_fd *fd, cram_container *c,
					   cram_block_compression_hdr *h);

/*! INTERNAL:
 * Encodes a slice compression header. 
 *
 * @return
 * Returns cram_block on success;
 *         NULL on failure
 */
cram_block *cram_encode_slice_header(cram_fd *fd, cram_slice *s);

/*! INTERNAL:
 * Encodes all slices in a container into blocks.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 *
 * FIXME: separate into encode_container and write_container. Ideally
 * we should be able to do read_container / write_container or
 * decode_container / encode_container.
 */
int cram_encode_container(cram_fd *fd, cram_container *c);

#ifdef __cplusplus
}
#endif

#endif
