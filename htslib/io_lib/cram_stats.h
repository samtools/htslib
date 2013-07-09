/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

#ifndef _CRAM_STATS_H_
#define _CRAM_STATS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "io_lib/hash_table.h"

cram_stats *cram_stats_create(void);
void cram_stats_add(cram_stats *st, int32_t val);
void cram_stats_del(cram_stats *st, int32_t val);
void cram_stats_dump(cram_stats *st);
void cram_stats_free(cram_stats *st);

/*
 * Computes entropy from integer frequencies for various encoding methods and
 * picks the best encoding.
 *
 * FIXME: we could reuse some of the code here for the actual encoding
 * parameters too. Eg the best 'k' for SUBEXP or the code lengths for huffman.
 *
 * Returns the best codec to use.
 */
enum cram_encoding cram_stats_encoding(cram_fd *fd, cram_stats *st);

#ifdef __cplusplus
}
#endif

#endif
