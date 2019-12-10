# HTSlib 64 bit reference positions

HTSlib version 1.10 onwards internally use 64 bit reference positions.  This
is to support analysis of species like axolotl, tulip and marbled lungfish
which have, or are expected to have,  chromosomes longer than two gigabases.

# File format support

Currently 64 bit positions can only be stored in SAM and VCF format files.
Binary BAM, CRAM and BCF cannot be used due to limitations in the formats
themselves.  As SAM and VCF are text formats, they have no limit on the
size of numeric values. Note that while 64 bit positions are supported by
default for SAM, for VCF they must be enabled explicitly at compile time
by editing Makefile and adding -DVCF_ALLOW_INT64=1 to CFLAGS.

# Compatibility issues to check

Various data structure members, function parameters, and return values have
been expanded from 32 to 64 bits.  As a result, some changes may be needed to
code that uses the library, even if it does not support long references.

## Variadic functions taking format strings

The type of various structure members (e.g. `bam1_core_t::pos`) and return
values from some functions (e.g. `bam_cigar2rlen()`) have been changed to
`hts_pos_t`, which is a 64-bit signed integer.  Using these in 32-bit
code will generally work (as long as the stored positions are within range),
however care needs to be taken when these values are passed directly
to functions like `printf()` which take a variable-length argument list and
a format string.

Header file `htslib/hts.h` defines macro `PRIhts_pos` which can be
used in `printf()` format strings to get the correct format specifier for
an `hts_pos_t` value.  Code that needs to print positions should be
changed from:

```c
printf("Position is %d\n", bam->core.pos);
```

to:

```c
printf("Position is %"PRIhts_pos"\n", bam->core.pos);
```

If for some reason compatibility with older versions of HTSlib (which do
not have `hts_pos_t` or `PRIhts_pos`) is needed, the value can be cast to
`int64_t` and printed as an explicitly 64-bit value:

```c
#include <inttypes.h> // For PRId64 and int64_t

printf("Position is %" PRId64 "\n", (int64_t) bam->core.pos);
```

Passing incorrect types to variadic functions like `printf()` can lead
to incorrect behaviour and security risks, so it important to track down
and fix all of the places where this may happen.  Modern C compilers like
gcc (version 3.0 onwards) and clang can check `printf()` and `scanf()`
parameter types for compatibility against the format string.  To
enable this, build code with `-Wall` or `-Wformat` and fix all the
reported warnings.

Where functions that take `printf`-style format strings are implemented,
they should use the appropriate gcc attributes to enable format string
checking.  `htslib/hts_defs.h` includes macros `HTS_FORMAT` and
`HTS_PRINTF_FMT` which can be used to provide the attribute declaration
in a portable way.  For example, `test/sam.c` uses them for a function
that prints error messages:

```
void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) fail(const char *fmt, ...) { /* ... */ }
```

## Implicit type conversions

Conversion of signed `int` or `int32_t` to `hts_pos_t` will always work.

Conversion of `hts_pos_t` to `int` or `int32_t` will work as long as the value
converted is within the range that can be stored in the destination.

Code that casts unsigned `uint32_t` values to signed with the expectation
that the result may be negative will no longer work as `hts_pos_t` can store
values over UINT32_MAX.  Such code should be changed to use signed values.

Functions hts_parse_region() and hts_parse_reg64() return special value
`HTS_POS_MAX` for regions which extend to the end of the reference.
This value is slightly smaller than INT64_MAX, but should be larger than
any reference that is likely to be used.  When cast to `int32_t` the
result should be `INT32_MAX`.

# Upgrading code to work with 64 bit positions

Variables used to store reference positions should be changed to
type `hts_pos_t`.  Use `PRIhts_pos` in format strings when printing them.

When converting positions stored in strings, use `strtoll()` in place of
`atoi()` or `strtol()` (which produces a 32 bit value on 64-bit Windows and
all 32-bit platforms).

Programs which need to look up a reference sequence length from a `sam_hdr_t`
structure should use `sam_hdr_tid2len()` instead of the old
`sam_hdr_t::target_len` array (which is left as 32-bit for reasons of
compatibility).  `sam_hdr_tid2len()` returns `hts_pos_t`, so works correctly
for large references.

Various functions which take pointer arguments have new versions which
support `hts_pos_t *` arguments.  Code supporting 64-bit positions should
use the new versions.  These are:

Original function  | 64-bit version
------------------ | --------------------
fai_fetch()        | fai_fetch64()
fai_fetchqual()    | fai_fetchqual64()
faidx_fetch_seq()  | faidx_fetch_seq64()
faidx_fetch_qual() | faidx_fetch_qual64()
hts_parse_reg()    | hts_parse_reg64() or hts_parse_region()
bam_plp_auto()     | bam_plp64_auto()
bam_plp_next()     | bam_plp64_next()
bam_mplp_auto()    | bam_mplp64_auto()

Limited support has been added for 64-bit INFO values in VCF files, for large
values in structural variant END tags.  New functions `bcf_update_info_int64()`
and `bcf_get_info_int64()` can be used to set and fetch 64-bit INFO values.
They both take arrays of `int64_t`.  `bcf_int64_missing` and
`bcf_int64_vector_end` can be used to set missing and vector end values in
these arrays.  The INFO data is stored in the minimum size needed, so there
is no harm in using these functions to store smaller integer values.

# Structure members that have changed size

```
File htslib/hts.h:
   hts_pair32_t::begin
   hts_pair32_t::end

   (typedef hts_pair_pos_t is provided as a better-named replacement for hts_pair32_t)

   hts_reglist_t::min_beg
   hts_reglist_t::max_end

   hts_itr_t::beg
   hts_itr_t::end
   hts_itr_t::curr_beg
   hts_itr_t::curr_end

File htslib/regidx.h:
   reg_t::start
   reg_t::end

File htslib/sam.h:
   bam1_core_t::pos
   bam1_core_t::mpos
   bam1_core_t::isize

File htslib/synced_bcf_reader.h:
   bcf_sr_regions_t::start
   bcf_sr_regions_t::end
   bcf_sr_regions_t::prev_start

File htslib/vcf.h:
   bcf_idinfo_t::info

   bcf_info_t::v1::i

   bcf1_t::pos
   bcf1_t::rlen
```

# Functions where parameters or the return value have changed size

Functions are annotated as follows:

* `[new]`  The function has been added since version 1.9
* `[parameters]` Function parameters have changed size
* `[return]` Function return value has changed size

```
File htslib/faidx.h:

   [new]        fai_fetch64()
   [new]        fai_fetchqual64()
   [new]        faidx_fetch_seq64()
   [new]        faidx_fetch_qual64()
   [new]        fai_parse_region()

File htslib/hts.h:

   [parameters] hts_idx_push()
   [new]        hts_parse_reg64()
   [parameters] hts_itr_query()
   [parameters] hts_reg2bin()

File htslib/kstring.h:

   [new]        kputll()

File htslib/regidx.h:

   [parameters] regidx_overlap()

File htslib/sam.h:

   [new]        sam_hdr_tid2len()
   [return]     bam_cigar2qlen()
   [return]     bam_cigar2rlen()
   [return]     bam_endpos()
   [parameters] bam_itr_queryi()
   [parameters] sam_itr_queryi()
   [new]        bam_plp64_next()
   [new]        bam_plp64_auto()
   [new]        bam_mplp64_auto()
   [parameters] sam_cap_mapq()
   [parameters] sam_prob_realn()

File htslib/synced_bcf_reader.h:

   [parameters] bcf_sr_seek()
   [parameters] bcf_sr_regions_overlap()

File htslib/tbx.h:

   [parameters] tbx_readrec()

File htslib/vcf.h:

   [parameters] bcf_readrec()
   [new]        bcf_update_info_int64()
   [new]        bcf_get_info_int64()
   [return]     bcf_dec_int1()
   [return]     bcf_dec_typed_int1()

```
