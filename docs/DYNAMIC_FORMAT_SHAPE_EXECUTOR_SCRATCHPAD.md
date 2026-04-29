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

## Open Questions

- How much of the gap is parse-loop dispatch versus generic encode cost?
- Can string width measurement be cached per shape region, or does row-local
  width variation force a cheap scan every time?
- Is it better to build several executor families by op sequence, or one generic
  fixed-shape executor with parse-function pointers?
- Do temporary fallback reason counters pay for themselves during iteration, or
  should they stay under an explicit debug environment variable?
