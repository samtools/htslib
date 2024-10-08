/*
Copyright (c) 2012-2013, 2018, 2024 Genome Research Ltd.
Author: James Bonfield <jkb@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! \file
 * Include cram.h instead.
 *
 * This is an internal part of the CRAM system and is automatically included
 * when you #include cram.h.
 *
 * Implements the decoding portion of CRAM I/O. Also see
 * cram_codecs.[ch] for the actual encoding functions themselves.
 */

#ifndef CRAM_DECODE_H
#define CRAM_DECODE_H

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
 * Loads and decodes the next slice worth of data.
 *
 * @return
 * Returns cram slice pointer on success;
 *         NULL on failure
 */
cram_slice *cram_next_slice(cram_fd *fd, cram_container **cp);

/*! INTERNAL:
 * Decode an entire slice from container blocks. Fills out s->crecs[] array.
 *
 * @return
 * Returns 0 on success;
 *        -1 on failure
 */
int cram_decode_slice(cram_fd *fd, cram_container *c, cram_slice *s,
                      sam_hdr_t *hdr);


/*! INTERNAL:
 * Converts a cram in-memory record into a bam in-memory record. We
 * pass a pointer to a bam_seq_t pointer along with the a pointer to
 * the allocated size. These can initially be pointers to NULL and zero.
 *
 * This function will reallocate the bam buffer as required and update
 * (*bam)->alloc accordingly, allowing it to be used within a loop
 * efficiently without needing to allocate new bam objects over and
 * over again.
 *
 * Returns the used size of the bam record on success
 *         -1 on failure.
 */
int cram_to_bam(sam_hdr_t *sh, cram_fd *fd, cram_slice *s,
                cram_record *cr, int rec, bam_seq_t **bam);

/*
 * Drains and frees the decode read-queue for a multi-threaded reader.
 */
void cram_drain_rqueue(cram_fd *fd);

#ifdef __cplusplus
}
#endif

#endif
