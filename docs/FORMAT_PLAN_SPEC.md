# FORMAT Plan Parser Spec

This document describes the intended direction for the experimental
`HTS_VCF_FORMAT_PLAN=1` VCF FORMAT parser.

## Goal

Keep the existing parser as the source of truth, but add a runtime-compiled
fast path for common FORMAT layouts.  The fast path should be opportunistic:
compile a plan for repeated FORMAT strings, execute known-safe operations
directly, and fall back to the generic parser whenever the record leaves the
supported subset.

## Architecture

The parser is tiered:

1. Exact kernels for dominant production layouts.  The current CCDG kernels cover
   `GT:AB:AD:DP:GQ:PL`, `GT:AD:DP:GQ:PL`,
   `GT:AB:AD:DP:GQ:PGT:PID:PL`, and `GT:AD:DP:GQ:PGT:PID:PL`.
2. A compiled op-list interpreter for regular FORMAT layouts.  It caches the
   FORMAT string, resolves header IDs once, then executes per-field operations
   for GT, integer vectors, float vectors, and strings.
3. Generic htslib parsing for everything else, including sample subsetting,
   duplicate FORMAT tags, undefined tags that require dummy header insertion,
   unsupported header types, malformed values, or future VCF constructs.

The cache key is the literal FORMAT column.  Record-specific widths are still
computed per row because BCF stores each FORMAT field as a rectangular
sample-by-value array, and the width depends on observed ploidy, vector length,
string length, and allele count.

## Correctness Rules

The planned parser must produce byte-identical BCF to the generic parser for any
record it claims.  If it cannot prove that, it must return `-3` so the existing
parser handles the record.

Required invariants:

- No planned parsing while `h->keep_samples` is active.
- Header IDs and types are resolved before execution.
- Duplicate tags use the generic parser.
- Undefined tags use the generic parser, preserving current dummy-header
  behavior and warnings.
- GT encoding must match generic htslib phasing semantics, including haploid
  genotypes and VCF 4.4 prefix phasing.
- Numeric vectors use observed row width and pad shorter samples with vector-end
  sentinels.
- Strings use observed maximum byte length and zero-pad shorter samples.
- Integer and float overflow/error behavior should either match generic htslib
  or force fallback.

## Current MVP

The current implementation keeps the CCDG exact kernels as the first tier and
adds a general compiled op-list tier for defined FORMAT fields with type
`String`, `Integer`, or `Float`.  The op-list tier handles:

- arbitrary field order,
- haploid, diploid, multidigit, missing, and phased GT values,
- integer and float vectors with row-local observed widths,
- string fields with row-local observed widths,
- multiallelic `Number=R` and `Number=G` rows by using observed vector width.

The MVP intentionally falls back for sample subsetting, duplicate tags,
undefined tags, unsupported header types, and malformed values.

After the row width pass, the interpreter resolves each cached FORMAT op to a
row-specific opcode such as `GT2`, `GT`, `INT1`, `INT2`, `INT3`, `INTN`,
`FLOAT1`, `FLOATN`, or `STR`.  This keeps layout coverage flexible while
memoizing the common "muscle memory" for repeated shapes.

For rows whose widths can be predicted from the header and allele count, the
interpreter can try a strict path before the observed-width pass.  The strict
path validates the observed maximum width while parsing and falls back to the
observed-width interpreter if the row is sparse, malformed, string-bearing, or
otherwise not byte-identical.  Common numeric opcode tapes such as
`GT2:INT2:INT1:INT1:INT3` and `GT2:FLOAT1:INT2:INT1:INT1:INT3` have
shape-level executors that avoid the per-op switch.

Planned integer parsing must be overflow-safe.  If a value is outside the BCF
int32 payload range, the planned parser falls back so the generic parser keeps
its warning and missing-value behavior.

## Edge Fixture

`test/format-plan-edge.vcf` is CCDG-shaped but includes records that exercise
common awkward cases:

- the exact CCDG layouts,
- reordered fields,
- multiallelic AD/PL and GL,
- haploid GT,
- multidigit allele indexes,
- fixed integer vectors,
- string FORMAT fields,
- exact-kernel fallbacks such as haploid GT and multidigit allele indexes.

Run:

```sh
./test/test_format_plan.sh
```

The script writes BCF through the generic parser and through
`HTS_VCF_FORMAT_PLAN=1`, compares them with `cmp`, and prints plan hit/fallback
statistics.  `HTS_VCF_FORMAT_PLAN=interp` or `HTS_VCF_FORMAT_PLAN=general`
skips the exact kernels and runs the compiled op-list interpreter directly,
which is useful for isolating interpreter performance.

## Next Work

- Add more exact kernels only after coverage data shows that they dominate real
  inputs.
- Add plan- or shape-level executors for dominant opcode sequences so hot rows
  can also avoid the per-sample opcode switch.
- Explore direct final-buffer output for validated fixed-width fields; this is
  likely higher leverage than adding more switch-level shape executors.
- Add overflow-compatible numeric parsing or force fallback before committing to
  the plan on extreme integer/float values.
- Integrate the edge fixture into the standard htslib test runner once the
  experimental flag graduates beyond local benchmarking.
