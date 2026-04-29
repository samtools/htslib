# Dynamic FORMAT Shape Executor Scratchpad

Date: 2026-04-29

Branch: `codex/vcf-avx-sanity`

## Goal

Make the general-purpose VCF FORMAT planned parser approach the handwritten
exact CCDG kernel speed without matching on field names.  The planner should stay
general, but the hot executor should become shape-specialized once a repeated
FORMAT layout proves stable.

The production htslib parser remains the source of truth.  Any optimized path
must either emit byte-identical BCF or return `-3` and fall back.

## Current Baseline

Known modes:

```sh
HTS_VCF_FORMAT_PLAN=0       # existing generic parser
HTS_VCF_FORMAT_PLAN=1       # exact CCDG kernels, then dynamic fallback
HTS_VCF_FORMAT_PLAN=interp  # dynamic planner only
HTS_VCF_FORMAT_PLAN_STATS=1 # counters from test/test_view
```

Current 10k CCDG sanity timing:

| Mode | VCF.gz read-only | VCF.gz -> uncompressed BCF |
|---|---:|---:|
| Baseline | 2.58 s | 2.83 s |
| Exact + dynamic fallback | 1.61 s | 1.86 s |
| Dynamic general only | 2.34 s | 2.55 s |

The performance target is the exact CCDG tier.  The first milestone is not to
delete exact kernels, but to make a dynamic shape executor selected without tag
name special cases reach the same neighborhood.

## Working Hypothesis

The exact kernels are faster because they do less work in the sample loop:

- no per-field switch dispatch for every sample,
- fewer scratch-buffer passes,
- direct writes into final BCF payloads for cheap fixed fields,
- integer min/max tracking is carried directly into BCF integer width selection,
- CCDG `PGT/PID` string handling is tailored instead of fully generic,
- fallback checks are simple and close to the parse step.

The dynamic planner should keep general discovery, but execute as a compact
fixed-shape kernel after resolving header metadata and row-local widths.

## Design Direction

Compile the literal FORMAT column plus header metadata into a plan as today, then
derive a shape descriptor from the row and header:

```text
GT2, FLOAT1, INT_R, INT1, INT1, STR1, STR1, INT_G
```

This shape says what to parse, not which tag names are present.  Field IDs,
header types, and BCF keys still come from the generic plan.

The executor should be monomorphic for common shapes:

- `GT2 + fixed numeric fields`
- `GT2 + FLOAT1 + fixed numeric fields`
- `GT2 + fixed numeric fields + fixed strings`
- `GT dynamic + fixed numeric fields`
- measured-width fallback for strings or sparse/non-fixed rows

The important constraint is to move per-field dispatch out of the per-sample hot
loop wherever possible.

## Correctness Rules

- Do not run planned parsing when `h->keep_samples` is active.
- Fall back on duplicate FORMAT tags, undefined tags, unsupported header types,
  malformed rows, unsupported GT shape, or integer/float behavior that cannot be
  made byte-identical.
- Preserve htslib GT semantics: haploid GT, missing alleles, multidigit allele
  indexes, phased/unphased state, and VCF 4.4 prefix phasing.
- Preserve vector-end padding and string zero-padding.
- Save and roll back `v->indiv.l` before any direct final-buffer write that may
  fall back.
- Keep exact kernels available as an oracle until dynamic shape execution closes
  the gap.

## Implementation Plan

1. Add shape classification to the dynamic general plan path.
   - Use existing `vcf_format_row_op_t` data where possible.
   - Recognize fixed-width rows derived from `Number=1`, `Number=R`,
     `Number=G`, and fixed `Number=N`.
   - Reject rows needing measured widths unless handled by a specific executor.

2. Add a first generic fixed-shape executor for CCDG-equivalent structures.
   - No tag-name matching.
   - Require leading `GT2`.
   - Support any mix/order of fixed INT/REAL fields with widths 1, R, G, or
     small fixed N.
   - Initially support `Number=1` strings with measured max width so `PGT/PID`
     can stay on the planned path.

3. Reduce hot-loop dispatch.
   - Precompute field offsets, widths, sizes, and parse actions.
   - Prefer executor-family loops over `switch (op->kind)` per field per sample.
   - Specialize common width parse helpers for 1, 2, 3, and small fixed widths.

4. Direct-write final BCF output when safe.
   - Continue direct `GT2` int8 output.
   - Track integer ranges while parsing and use known-range encoding.
   - For floats, serialize directly when field width is fixed.
   - For strings, write final char payload after width is known.

5. Instrument fallback reasons and executor choices.
   - Add temporary counters or debug logging gated by env vars if useful.
   - Track shape hits, shape fallbacks, strict hits, measured fallback hits.

6. Benchmark and iterate.
   - Correctness: `./test/test_format_plan.sh`
   - CCDG subset: `/tmp/ccdg_chr22_10k.vcf.gz`
   - Compare baseline, exact, and dynamic-only output with `cmp`.
   - Primary target: dynamic-only `HTS_VCF_FORMAT_PLAN=interp` approaching
     exact-mode time on VCF.gz -> uncompressed BCF.

## Test Data

Full source:

```text
/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/data/original/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

Benchmark subset:

```text
/tmp/ccdg_chr22_10k.vcf.gz
/tmp/ccdg_chr22_10k.bcf
```

Correctness fixture:

```text
test/format-plan-edge.vcf
```

## Current Scratch Notes

- `HTS_VCF_FORMAT_PLAN=interp` is the key mode for dynamic executor progress.
- Exact kernels should remain until dynamic-only is close enough to make them
  redundant.
- Avoid hardcoding `AD`, `PL`, `DP`, `GQ`, `AB`, `PGT`, or `PID`; use their
  header-derived type/number/width instead.
- CCDG-like FORMAT distributions are still the first target because they provide
  a real, repeatable workload and a clear oracle.

## 2026-04-29 Iteration Notes

Implemented the first dynamic likelihood-shape executor in `vcf.c`.

What changed:

- Added an optional `HTS_VCF_FORMAT_PLAN_SHAPE_STATS` counter path in
  `test/test_view`.
- Relaxed strict string handling so `Type=String,Number=1` FORMAT fields can be
  handled by planned parsing with row-local byte-width measurement.
- Added a shape-specific width derivation for CCDG-like layouts where `AD` may
  be declared as `Number=.` in the header but the row shape proves the observed
  width is `n_allele`.
- Added a straight-line dynamic executor for:

```text
GT2, optional FLOAT1, INT[n_allele], INT1, INT1,
optional STR1, optional STR1, INT[n_allele * (n_allele + 1) / 2]
```

This executor is selected by FORMAT type/order/width, not by tag names.  It
still validates observed AD/PL counts and falls back on mismatch.

Latest 10k CCDG VCF.gz -> uncompressed BCF single-run timings on the rebuilt
worktree:

| Mode | Wall | User | Notes |
|---|---:|---:|---|
| Baseline | 2.78 s | 2.56 s | `HTS_VCF_FORMAT_PLAN=0` |
| Exact CCDG | 1.78 s | 1.61 s | exact kernels, shape hits 0 |
| Dynamic shape | 2.53 s | 1.71 s | `interp`, shape hits 10,000 |

`cmp` passed for both dynamic-shape and exact outputs against baseline BCF.

The important result is CPU parity is close: dynamic shape is within about 6% of
exact on user time in this run.  Wall time is noisier, likely output/cache
effects, and should be rerun in a tighter benchmark loop.

Next likely cuts:

- Cache shape classification per `(header, FORMAT)` plan so we do less
  per-record type/order checking.
- Split phase and non-phase shape executors to remove `has_phase` branches from
  the sample loop.
- Consider separate `has_float` executor variants for the same reason.
- Compare a no-shape-stats build/run to estimate counter overhead, though it is
  probably minor.
- Once dynamic shape is consistently at parity, demote the exact CCDG kernels to
  oracle-only or remove them.

## 2026-04-29 Large-Corpus Check

The small public/synthetic slices were too short to provide timing signal, so
the meaningful benchmark set moved to `bench/format-shape/large/inputs.tsv`.
The large corpus includes the CCDG 10k subset, full 1000 Genomes chr22
genotypes, and 2,048-sample generated FORMAT workloads.

Latest large run used:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All exact/interp outputs compared byte-identical to baseline.  Timing summary:

| Input | Baseline user | Exact user | Dynamic interp user | Shape hits |
|---|---:|---:|---:|---:|
| CCDG 10k | 2.60 s | 1.64 s | 1.63 s | 10,000 |
| 1000G chr22 full GT | 26.76 s | 9.32 s | 9.32 s | 0 |
| Large CCDG-like synthetic | 4.10 s | 2.77 s | 2.87 s | 20,000 |
| Large reordered likelihood | 2.97 s | 2.65 s | 2.62 s | 0 |
| Large multiallelic likelihood | 3.22 s | 2.09 s | 2.05 s | 16,000 |
| Large float/string | 2.95 s | 2.84 s | 2.84 s | 0 |

The dynamic likelihood shape path is now at parity or close enough on the
meaningful workloads.  The remaining visible gap is the generated CCDG-like
phase-heavy synthetic case, where dynamic-only is about 3-4% slower than exact.
That looks acceptable for this checkpoint; the next optimization target remains
cached shape classification to remove repeated deterministic row-level checks.

## 2026-04-29 Cached Shape Classification

Added FORMAT-level likelihood-shape classification to the dynamic general plan.
The cache only records deterministic facts from the FORMAT/header order and
types:

```text
GT2, optional FLOAT1, INT[n_allele], INT1, INT1,
optional STR1, optional STR1, INT[ploidy likelihood width]
```

Row-level facts remain uncached.  Each record still validates `n_allele`,
AD/PL widths, GT syntax, observed vector counts, separators, sample count, and
phase-string widths before using the likelihood executor.

The large benchmark corpus now includes four extra cache-regression workloads:

- variable-width `PGT/PID` likelihood rows,
- likelihood rows with mixed row-local fallbacks and later positive hits,
- GT-first but wrong-order likelihood-like rows,
- non-likelihood rows with two strings plus float vectors.

Latest full large-corpus run:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All exact and interp outputs compared byte-identical to baseline.  Highlights:

| Input | Exact user | Dynamic interp user | Dynamic shape attempts | Dynamic shape hits |
|---|---:|---:|---:|---:|
| CCDG 10k | 1.61 s | 1.60 s | 10,000 | 10,000 |
| 1000G chr22 full GT | 9.16 s | 9.11 s | 0 | 0 |
| Large CCDG-like synthetic | 2.74 s | 2.69 s | 20,000 | 20,000 |
| Large multiallelic likelihood | 2.05 s | 1.99 s | 16,000 | 16,000 |
| Variable phase widths | 2.00 s | 1.99 s | 12,000 | 12,000 |
| Mixed row-local fallbacks | 1.56 s | 1.58 s | 11,295 | 10,236 |
| GT-first reordered negative | 1.50 s | 1.50 s | 0 | 0 |
| Two-string float negative | 2.36 s | 2.32 s | 0 | 0 |

The important negative-cache result is the full 1000G GT-only workload:
dynamic mode no longer pays 1,103,547 failed likelihood-shape probes.

## 2026-04-29 GT-Only Fast Path

Added a tiny general-plan executor for the common `FORMAT=GT` / diploid `GT2`
shape.  This is still shape-based rather than data-set specific:

- requires a single FORMAT op and that op must be `GT`,
- requires allele indexes that fit the existing one-digit `GT2` parser,
- falls through to the existing strict/measured paths for haploid, dynamic GT,
  malformed rows, or any unsupported row-local detail.

Also tightened the older exact-name CCDG kernels so they only claim a FORMAT
after checking the relevant header types and scalar counts.  A new
`format-plan-header-mismatch.vcf` fixture keeps this honest by using CCDG-shaped
names with `AD` declared as a string.

Latest full large-corpus run remained byte-identical to baseline for exact and
interp.  The main win is the full 1000G chr22 GT-only workload:

| Input | Baseline user | Exact user | Dynamic interp user |
|---|---:|---:|---:|
| 1000G chr22 full GT | 24.86 s | 5.77 s | 5.61 s |

The previous cached-shape run was about 9.1 s user in dynamic mode on this
input, so the direct GT-only executor removes roughly 39% of the remaining
planned-parser CPU for this large real workload.

## 2026-04-29 Multiallelic Parse Tightening

Added two small low-risk likelihood executor refinements:

- avoid retrying the likelihood shape executor inside the strict path when the
  same row already reached the likelihood executor and failed row-local checks;
- add fixed-width integer vector parsers for AD width 4 and PL widths 6/10,
  covering common triallelic and quad-allelic `Number=G` likelihood rows.

The fixed-width parsers still use the same scalar integer parser and range
tracking, and they preserve short-vector padding and trailing-comma fallback
behavior.

Small edge coverage now includes:

- row-local likelihood fallback from short AD/PL in individual samples,
- missing AD/PL with another sample proving full row width,
- GT-only fast-path hits plus haploid and multidigit GT fallbacks.

Latest full large-corpus run remained byte-identical to baseline.  Timings were
noisier overall than the previous pass, but the important rows were:

| Input | Exact user | Dynamic interp user | Notes |
|---|---:|---:|---|
| CCDG 10k | 1.59 s | 1.56 s | likelihood shape parity |
| 1000G chr22 full GT | 5.64 s | 5.68 s | GT-only fast path retained |
| Large multiallelic likelihood | 2.26 s | 2.07 s | dynamic ahead of exact |
| Mixed row-local fallbacks | 1.72 s | 1.74 s | byte-clean fallback path |

## 2026-04-29 No-Special Integer Encode

Added a conservative `has_special` bit to planned integer range tracking.  The
parser now records when it has observed `bcf_int32_missing` or
`bcf_int32_vector_end`, including vector-end padding from short fixed-width
vectors.  The known-range encoder uses that proof to skip per-value sentinel
checks in int8/int16 output only when the row contains no missing/vector-end
values.

Safety rule: min/max alone never proves this.  Missing and vector-end sentinels
can still select int8/int16 BCF encodings, so the fast loop is gated only by the
parser-maintained flag.

Small edge coverage now includes integer boundary rows spanning int8/int16/int32
choices, plus existing rows with scalar missing values, short fixed vectors, and
explicit vector missing values.

Latest full large-corpus run:

| Input | Exact user | Dynamic interp user | Notes |
|---|---:|---:|---|
| CCDG 10k | 1.74 s | 1.73 s | real likelihood parity |
| 1000G chr22 full GT | 6.02 s | 6.09 s | GT-only path retained |
| Large CCDG-like synthetic | 3.03 s | 2.98 s | dynamic slightly ahead |
| Large multiallelic likelihood | 2.29 s | 2.12 s | dynamic ahead |
| Mixed row-local fallbacks | 1.71 s | 1.72 s | byte-clean fallback path |

All exact and interp outputs compared byte-identical to baseline.

## Open Questions

- How much of the gap is parse-loop dispatch versus generic encode cost?
- Can string width measurement be cached per shape region, or does row-local
  width variation force a cheap scan every time?
- Is it better to build several executor families by op sequence, or one generic
  fixed-shape executor with parse-function pointers?
- Do temporary fallback reason counters pay for themselves during iteration, or
  should they stay under an explicit debug environment variable?
