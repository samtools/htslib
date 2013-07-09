/*! \file
 * CRAM interface.
 *
 * Consider using the higher level scram_*() API for programs that wish to
 * be file format agnostic.
 *
 * This API should be used for CRAM specific code. The specifics of the
 * public API are implemented in cram_io.h, cram_encode.h and cram_decode.h
 * although these should not be included directly (use this file instead).
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef SAMTOOLS
#  include "io_lib/cram_samtools.h"
#endif

#ifndef _CRAM_H_
#define _CRAM_H_

#include "io_lib/sam_header.h"
#include "cram_structs.h"
#include "cram_io.h"
#include "cram_encode.h"
#include "cram_decode.h"
#include "cram_stats.h"
#include "cram_codecs.h"
#include "cram_index.h"

#endif

#ifdef __cplusplus
}
#endif
