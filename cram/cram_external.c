/*
Copyright (c) 2015 Genome Research Ltd.
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
 * External CRAM interface.
 *
 * Internally we're happy to use macros and to grub around in the cram
 * structures.  This isn't very sustainable for an externally usable
 * ABI though, so we have anonymous structs and accessor functions too
 * to permit software such as samtools reheader to manipulate cram
 * containers and blocks in a robust manner.
 */

#include "htslib/hfile.h"
#include "cram/cram.h"

/*
 *-----------------------------------------------------------------------------
 * cram_fd
 */
SAM_hdr *cram_fd_get_header(cram_fd *fd) { return fd->header; }
void cram_fd_set_header(cram_fd *fd, SAM_hdr *hdr) { fd->header = hdr; }

int cram_fd_get_version(cram_fd *fd) { return fd->version; }
void cram_fd_set_version(cram_fd *fd, int vers) { fd->version = vers; }

int cram_major_vers(cram_fd *fd) { return CRAM_MAJOR_VERS(fd->version); }
int cram_mINOR_vers(cram_fd *fd) { return CRAM_MINOR_VERS(fd->version); }

hFILE *cram_fd_get_fp(cram_fd *fd) { return fd->fp; }
void cram_fd_set_fp(cram_fd *fd, hFILE *fp) { fd->fp = fp; }


/*
 *-----------------------------------------------------------------------------
 * cram_container
 */
int32_t cram_container_get_length(cram_container *c) {
    return c->length;
}

void cram_container_set_length(cram_container *c, int32_t length) {
    c->length = length;
}


int32_t cram_container_get_num_blocks(cram_container *c) {
    return c->num_blocks;
}

void cram_container_set_num_blocks(cram_container *c, int32_t num_blocks) {
    c->num_blocks = num_blocks;
}


/* Returns the landmarks[] array and the number of elements
 * in num_landmarks.
 */
int32_t *cram_container_get_landmarks(cram_container *c, int32_t *num_landmarks) {
    *num_landmarks = c->num_landmarks;
    return c->landmark;
}

/* Sets the landmarks[] array (pointer copy, not a memory dup) and
 * num_landmarks value.
 */
void cram_container_set_landmarks(cram_container *c, int32_t num_landmarks,
				  int32_t *landmarks) {
    c->num_landmarks = num_landmarks;
    c->landmark = landmarks;
}


/*
 *-----------------------------------------------------------------------------
 * cram_block
 */
int32_t cram_block_get_content_id(cram_block *b)  { return b->content_id; }
int32_t cram_block_get_comp_size(cram_block *b)   { return b->comp_size; }
int32_t cram_block_get_uncomp_size(cram_block *b) { return b->uncomp_size; }
int32_t cram_block_get_crc32(cram_block *b)       { return b->crc32; }
void *  cram_block_get_data(cram_block *b)        { return b->data; }

void cram_block_set_content_id(cram_block *b, int32_t id) { b->content_id = id; }
void cram_block_set_comp_size(cram_block *b, int32_t size) { b->comp_size = size; }
void cram_block_set_uncomp_size(cram_block *b, int32_t size) { b->uncomp_size = size; }
void cram_block_set_crc32(cram_block *b, int32_t crc) { b->crc32 = crc; }
void cram_block_set_data(cram_block *b, void *data) { b->data = data; }

int cram_block_append(cram_block *b, void *data, int size) {
    BLOCK_APPEND(b, data, size);
    return BLOCK_DATA(b) ? 0 : -1; // It'll do for now...
}
void cram_block_update_size(cram_block *b) { BLOCK_UPLEN(b); }

// Offset is known as "size" internally, but it can be confusing.
size_t cram_block_get_offset(cram_block *b) { return BLOCK_SIZE(b); }
void cram_block_set_offset(cram_block *b, size_t offset) { BLOCK_SIZE(b) = offset; }
