# Dynamic FORMAT Plan: Current Implementation

This document describes the implementation currently present in `vcf.c`, the
correctness boundaries, and the latest benchmark results.

## Entry Point

`vcf_parse_format()` first calls `vcf_parse_format_planned()` when
`HTS_VCF_FORMAT_PLAN` is enabled.  The planned path either parses the whole
FORMAT column or returns `-3`, allowing the production parser to handle the
column unchanged.

```text
HTS_VCF_FORMAT_PLAN enabled
  -> compile FORMAT/header to per-tag plan
  -> resolve row-local widths
  -> composable executor
  -> production fallback on unsupported or suspicious rows
```

Enabled spellings are `1`, `interp`, and `general`; all route through the same
dynamic executor.

## Plan Compilation

Plans are cached by header pointer plus literal FORMAT string.  This is
important because VCF header IDs, declared types, and number models are
header-local.

The compile step rejects:

- undefined FORMAT tags;
- duplicate FORMAT tags;
- unsupported header types;
- unsupported number models;
- `GT` declarations that are not `Type=String,Number=1`.

Undefined tags intentionally fall back to the production parser so existing
dummy-header repair and warning behavior is preserved.

## Supported Operations

The current executor supports:

- `GT`, with fast `GT2` storage when the row is diploid and simple;
- integer fields with fixed `Number=N`, `Number=A`, `Number=R`, `Number=G`, or
  bounded measured `Number=.` row widths;
- float fields with the same number models as integer fields;
- string fields declared as `Type=String,Number=1`, measured per row.

Header-derived widths are resolved per record.  `Number=A`, `Number=R`, and
`Number=G` depend on the current allele count.  String and `Number=.` numeric
fields use a row-local measurement pass.

## Executor

BCF stores FORMAT data transposed by tag: all samples for FORMAT op 0, then all
samples for FORMAT op 1, and so on.  The dynamic executor parses VCF samples in
sample-major order and writes that transposed BCF layout.

Leading fixed-width `GT2` and `FLOAT1` rows can write directly into `v->indiv`.
Other rows are staged in header scratch memory, then encoded after sample
parsing so integer range and observed-width metadata are known.

For fixed-width vector fields, the executor can compact underfilled rows to the
observed row maximum before BCF encoding.  This avoids whole-row fallback when
the production parser would also emit a narrower byte-identical vector width.

## Guard Policy

Each cached dynamic plan has a small runtime guard:

- attempts, hits, fallbacks;
- consecutive miss streak;
- temporary cooldown.

An isolated fallback does not disable the fast path.  A plan is paused after
eight consecutive misses, or after at least 128 attempts with more than 10%
fallbacks.  After 256 skipped records, the plan probes again so later stable
regions can recover the optimized path.

## Correctness Rules

The planned parser must preserve these invariants:

- no planned parsing while `h->keep_samples` is active;
- header IDs, types, and number models are resolved before execution;
- duplicate or undefined tags use the production parser;
- unsupported GT encodings force fallback;
- numeric vectors preserve observed width and vector-end padding;
- strings use observed maximum byte length and zero-pad shorter samples;
- integer and float overflow/error behavior must match production htslib or
  force fallback;
- direct writes to `v->indiv` must roll back before fallback.

Focused validation lives in `./test/test_format_plan.sh`.  It compares
production parsing, `HTS_VCF_FORMAT_PLAN=1`, and the `interp` alias byte-for-byte
with `cmp`.

## Current Source Delta

After removing the old exact kernels and SIMD tab scanner, the live parser/test
hook delta relative to `origin/develop` is:

| File | Added lines |
|---|---:|
| `vcf.c` | 1,467 |
| `test/test_view.c` | 14 |

## Large Corpus Benchmark

Command:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-dynamic-trim-plan \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All planned outputs compared byte-identical to baseline.

| Input | Baseline user | Plan user | Hits/fallback |
|---|---:|---:|---:|
| CCDG 10k | 2.62 s | 2.25 s | 8,396 / 1,604 |
| 1000G chr22 full GT | 26.05 s | 7.98 s | 1,103,547 / 0 |
| Large CCDG-like synthetic | 4.24 s | 3.78 s | 20,000 / 0 |
| Large reordered likelihood | 3.00 s | 2.42 s | 20,000 / 0 |
| Large multiallelic likelihood | 3.16 s | 2.73 s | 16,000 / 0 |
| Large float/string | 2.93 s | 2.97 s | 16,000 / 0 |
| Variable phase widths | 2.61 s | 2.50 s | 12,000 / 0 |
| Mixed row-local fallbacks | 2.22 s | 1.87 s | 12,000 / 0 |
| GT-first reordered negative | 1.75 s | 1.44 s | 12,000 / 0 |
| Two-string float negative | 2.28 s | 2.56 s | 12,000 / 0 |

## bcftools Production-Style Benchmark

A clean bcftools `develop` worktree was built at:

```text
/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan
```

using this htslib worktree:

```sh
make HTSDIR=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/htslib-vcf-avx-sanity bcftools
```

Timing command:

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools \
  bench/format-shape/scripts/run_bcftools_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

The runner uses `bcftools view --no-version -Ob -l 0 [--threads N]`.  All
planned outputs compared byte-identical to baseline.

| Input | Threads | Baseline real | Plan real | Speedup | Baseline user | Plan user |
|---|---:|---:|---:|---:|---:|---:|
| 1000G chr22 full GT | 0 | 27.48 s | 8.99 s | 3.06x | 25.94 s | 8.05 s |
| 1000G chr22 full GT | 2 | 26.59 s | 6.99 s | 3.80x | 28.82 s | 9.04 s |
| 1000G chr22 full GT | 4 | 26.71 s | 6.94 s | 3.85x | 28.83 s | 9.08 s |
| 1000G chr22 full GT | 8 | 26.62 s | 6.96 s | 3.82x | 28.71 s | 9.38 s |
| Large CCDG-like synthetic | 0 | 4.43 s | 3.94 s | 1.12x | 4.11 s | 3.66 s |
| Large CCDG-like synthetic | 2 | 3.46 s | 3.01 s | 1.15x | 4.50 s | 4.06 s |
| Large CCDG-like synthetic | 4 | 3.47 s | 3.02 s | 1.15x | 4.51 s | 4.09 s |
| Large CCDG-like synthetic | 8 | 3.46 s | 3.00 s | 1.15x | 4.50 s | 4.05 s |

## Interpretation

The dynamic path gives a large production-visible win for sample-rich GT-only
VCFs.  On likelihood-heavy rows, it is consistently faster but still limited by
generic per-op work, string/width handling, and IO/compression costs.  Some
float/string-heavy layouts remain near parity or slightly slower than baseline.

## Remaining Work

- Add selected-sample support so `keep_samples` does not force fallback.
- Reduce per-sample opcode dispatch in hot FORMAT layouts.
- Improve string and measured-width handling without losing byte identity.
- Consider a later executor-generation layer if generic per-op dispatch remains
  the main gap to historical exact-kernel speed.
