# FORMAT Planner Overview

The FORMAT planner is an optional fast path for parsing VCF sample columns into
BCF.  It lives inside `vcf.c` and is disabled by default.  Set
`HTS_VCF_FORMAT_PLAN=1` to try the planned path before the normal generic
FORMAT parser.

The planner is deliberately conservative: if a FORMAT column is unsupported,
malformed, tied to a dirty header, or otherwise suspicious, it falls back for
the whole column and lets the existing parser own the result.  Byte-identical
BCF output versus the generic parser is the hard correctness requirement.

## Control Surface

- `HTS_VCF_FORMAT_PLAN=1` enables the planned FORMAT parser.
- Unset, `0`, and any other value use the generic parser only.
- `HTS_VCF_FORMAT_PLAN_STATS=1` prints diagnostic counters at process exit
  when the planner is enabled and reached.

The feature does not add public API.  Planner state is stored in private
`bcf_hdr_aux_t` data and is invalidated when `bcf_hdr_sync()` rebuilds header
dictionaries.

## Architecture

The planned path has four stages:

1. Cache lookup or compilation for the literal FORMAT string under the active
   header generation.
2. Row-local width resolution for `Number=A/R/G`, bounded `Number=.`, and
   string fields that require scanning sample text.
3. Sample-major parsing into transposed FORMAT rows.
4. BCF row encoding using the same public record layout as the generic parser.

Compiled plans describe tags at FORMAT-field granularity, not as exact whole
string kernels.  That lets the same executor handle layouts such as `GT:AD`,
`GT:AD:DP:PL`, reordered fields, and additional supported header-described
tags without adding one hand-written parser per FORMAT string.

## Supported Shapes

The planned executor currently covers:

- `GT` for the current `GT2` subset: two single-character alleles or missing
  values separated by `/` or `|`, such as `0/1`, `.|.`, `0|.`, and `.|0`.
- Integer scalar fields.
- Integer vector fields through a single shared `INTVEC` path.
- Float scalar and vector fields.
- Fixed `Number=1` strings.
- Header-derived `Number=A`, `Number=R`, and `Number=G` numeric widths, capped
  by the planner's numeric width limit.
- Bounded measured `Number=.` numeric rows and bounded measured strings.
- Selected-sample parsing via `bcf_hdr_set_samples()`.

The planner falls back on undefined tags, duplicate FORMAT tags, unsupported
header type/number models, dirty headers, unsupported GT shapes, malformed
separators, unsafe row widths, and unprofitable string/float-heavy layouts.
Haploid, polyploid, multi-digit-allele, and very high allele-count GT shapes
therefore stay on the generic parser.  Measured-string plus float-vector
compositions also fall back unless the FORMAT row contains integer-vector work
that makes the planned path worthwhile.

## Current Simplification

The integer-vector executor was simplified in
`codex/simplify-format-plan-executor`.  Width-specific row kinds and parsers
for `INT2`, `INT3`, and selected fixed widths were removed.  The executor now
uses:

- `GT2`
- `INT1`
- `INTVEC`
- `FLOAT1`
- `FLOATN`
- `STR`

This removed about 283 lines from `vcf.c` and made the code less shaped like a
collection of hand-specialized kernels.  Local untracked corpus-snapshot
benchmarking showed the simplified executor was slightly faster overall, with
byte-identical output against the pre-simplification planner.
