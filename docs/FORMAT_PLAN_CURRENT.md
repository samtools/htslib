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
  -> fetch or compile header-owned FORMAT/header plan
  -> resolve row-local widths
  -> composable executor
  -> production fallback on unsupported or suspicious rows
```

Enabled spellings are `1`, `interp`, and `general`; all route through the same
dynamic executor.

## Plan Compilation

Plans are cached in private `bcf_hdr_aux_t` state by literal FORMAT string plus
the header's private FORMAT-plan generation.  This is important because VCF
header IDs, declared types, and number models are header-local.  The cache grows
from 16 entries up to 128 entries, uses heap storage for long FORMAT strings,
and also caches unsupported schemas so repeated odd rows do not repeatedly pay
compile cost.

`bcf_hdr_sync()` clears the header-owned plan cache and increments the private
generation after header dictionaries are rebuilt.  The planner also refuses to
compile while `h->dirty` is set, leaving unsynced or header-repair cases on the
production parser.

The cache and per-plan guard counters are mutable header-owned state, like other
htslib header scratch storage.  Callers should not concurrently parse through
the same `bcf_hdr_t` from multiple threads.

The compile step rejects:

- undefined FORMAT tags;
- duplicate FORMAT tags;
- unsupported header types;
- unsupported number models;
- `GT` declarations that are not `Type=String,Number=1`.
- string-plus-float-vector schemas with too little integer-vector work to repay
  the dynamic path's width-measurement cost.

Undefined tags intentionally fall back to the production parser so existing
dummy-header repair and warning behavior is preserved.

## Supported Operations

The current executor supports:

- `GT`, with fast `GT2` storage when the row is diploid and simple;
- integer fields with fixed `Number=N`, `Number=A`, `Number=R`, `Number=G`, or
  bounded measured `Number=.` row widths;
- float fields with the same number models as integer fields;
- string fields declared as `Type=String,Number=1`, measured per row.
- `bcf_hdr_set_samples()` / `keep_samples`, by scanning the original sample
  columns and writing only retained samples densely into the planned BCF output.

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

- no planned parsing while the header has unsynced dictionary changes;
- header IDs, types, and number models are resolved before execution;
- selected-sample parsing must honor `h->keep_samples`, use `h->nsamples_ori`
  for input-column scans, and set `v->n_sample` to the retained sample count;
- duplicate or undefined tags use the production parser;
- unprofitable string/float-heavy schemas use the production parser;
- unsupported GT encodings force fallback;
- numeric vectors preserve observed width and vector-end padding;
- strings use observed maximum byte length and zero-pad shorter samples;
- integer and float overflow/error behavior must match production htslib or
  force fallback;
- direct writes to `v->indiv` must roll back before fallback.

Focused validation lives in `./test/test_format_plan.sh`.  It compares
production parsing, `HTS_VCF_FORMAT_PLAN=1`, and the `interp` alias byte-for-byte
with `cmp`.  The script also checks selected-sample parsing for explicit
inclusion and exclusion lists (`S1,S3`, `S2`, and `^S2`).  `test/format-plan-cache.vcf`
additionally exercises more than 16 distinct FORMAT schemas and a literal FORMAT
string longer than the old fixed cache key.  `test/test_format_plan_cache`
mutates and resyncs a header after a plan has been compiled for the same FORMAT
string, then verifies the row is planned again with the new metadata.

## Current Source Delta

After removing the old exact kernels and SIMD tab scanner, then hardening the
dynamic cache, the live parser/test hook delta relative to `origin/develop` is:

| File | Added lines |
|---|---:|
| `vcf.c` | 1,703 |
| `Makefile` | 6 |
| `test/test_format_plan.sh` | 48 |
| `test/test_format_plan_cache.c` | 130 |
| `test/test_view.c` | 23 |
| `test/format-plan-cache.vcf` | 61 |
| `test/format-plan-profitability.vcf` | 8 |

## Large Corpus Benchmark

Command:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-profit-gate \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All planned outputs compared byte-identical to baseline.

| Input | Baseline user | Plan user | Hits/fallback |
|---|---:|---:|---:|
| CCDG 10k | 2.47 s | 2.15 s | 8,396 / 1,604 |
| 1000G chr22 full GT | 25.25 s | 7.82 s | 1,103,547 / 0 |
| Large CCDG-like synthetic | 4.02 s | 3.64 s | 20,000 / 0 |
| Large reordered likelihood | 2.95 s | 2.40 s | 20,000 / 0 |
| Large multiallelic likelihood | 3.15 s | 2.76 s | 16,000 / 0 |
| Large float/string | 2.96 s | 2.89 s | 0 / 16,000 |
| Variable phase widths | 2.60 s | 2.46 s | 12,000 / 0 |
| Mixed row-local fallbacks | 2.19 s | 1.84 s | 12,000 / 0 |
| GT-first reordered negative | 1.72 s | 1.37 s | 12,000 / 0 |
| Two-string float negative | 2.29 s | 2.26 s | 0 / 12,000 |

## Full Threaded Corpus Benchmark

Command:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-threaded-profit-gate \
  bench/format-shape/scripts/run_thread_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

All 40 planned outputs compared byte-identical to baseline.  Detailed timings
are in `bench/format-shape/large/results-threaded-profit-gate/timings.tsv`; the
table below summarizes real-time speedup.

| Input | 0 threads | 2 threads | 4 threads | 8 threads |
|---|---:|---:|---:|---:|
| CCDG 10k | 1.13x | 1.15x | 1.16x | 1.15x |
| 1000G chr22 full GT | 3.10x | 3.73x | 4.34x | 3.88x |
| Large CCDG-like synthetic | 1.12x | 1.14x | 1.13x | 1.13x |
| Large reordered likelihood | 1.23x | 1.33x | 1.32x | 1.29x |
| Large multiallelic likelihood | 1.16x | 1.22x | 1.22x | 1.22x |
| Large float/string | 1.01x | 0.97x | 1.04x | 1.00x |
| Variable phase widths | 1.06x | 1.10x | 1.11x | 1.09x |
| Mixed row-local fallbacks | 1.18x | 1.25x | 1.31x | 1.23x |
| GT-first reordered negative | 1.22x | 1.31x | 1.32x | 1.32x |
| Two-string float negative | 1.00x | 1.00x | 1.01x | 1.00x |

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

`bench/format-shape/large/threaded-inputs.tsv` now mirrors the full large
corpus from `large/inputs.tsv`, so threaded runs cover all real and synthetic
workload shapes rather than only the earlier two representative rows.

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

## bcftools Selected-Sample Benchmark

The same bcftools runner can select the first N samples from each input with
`SAMPLE_COUNT=N`.  This exercises the `bcf_hdr_set_samples()` / `keep_samples`
path through bcftools rather than only through the test harness.

Command:

```sh
BCFTOOLS=/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/bcftools-htslib-vcf-plan/bcftools \
SAMPLE_COUNT=2 KEEP_OUTPUTS=0 \
OUTDIR=bench/format-shape/large/results-bcftools-keep2 \
  bench/format-shape/scripts/run_bcftools_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

All 40 planned outputs compared byte-identical to baseline.  The table shows
real-time and user-time speedup for selecting two samples from every input that
has samples; sites-only inputs naturally run without `-s`.

| Input | Threads | Real speedup | User speedup |
|---|---:|---:|---:|
| CCDG 10k | 0 | 1.12x | 1.12x |
| CCDG 10k | 2 | 1.12x | 1.11x |
| CCDG 10k | 4 | 1.13x | 1.12x |
| CCDG 10k | 8 | 1.11x | 1.10x |
| 1000G chr22 full GT | 0 | 2.71x | 2.73x |
| 1000G chr22 full GT | 2 | 2.83x | 2.44x |
| 1000G chr22 full GT | 4 | 2.94x | 2.52x |
| 1000G chr22 full GT | 8 | 3.06x | 2.61x |
| Large CCDG-like synthetic | 0 | 1.07x | 1.08x |
| Large CCDG-like synthetic | 2 | 1.10x | 1.07x |
| Large CCDG-like synthetic | 4 | 1.09x | 1.07x |
| Large CCDG-like synthetic | 8 | 1.09x | 1.07x |
| Large reordered likelihood | 0 | 1.15x | 1.17x |
| Large reordered likelihood | 2 | 1.22x | 1.15x |
| Large reordered likelihood | 4 | 1.23x | 1.17x |
| Large reordered likelihood | 8 | 1.22x | 1.16x |
| Large multiallelic likelihood | 0 | 1.13x | 1.13x |
| Large multiallelic likelihood | 2 | 1.14x | 1.11x |
| Large multiallelic likelihood | 4 | 1.16x | 1.12x |
| Large multiallelic likelihood | 8 | 1.18x | 1.13x |
| Large float/string | 0 | 1.02x | 1.01x |
| Large float/string | 2 | 0.99x | 0.99x |
| Large float/string | 4 | 1.01x | 1.00x |
| Large float/string | 8 | 0.97x | 0.98x |
| Variable phase widths | 0 | 1.04x | 1.05x |
| Variable phase widths | 2 | 1.05x | 1.05x |
| Variable phase widths | 4 | 1.05x | 1.04x |
| Variable phase widths | 8 | 1.06x | 1.05x |
| Mixed row-local fallbacks | 0 | 1.14x | 1.16x |
| Mixed row-local fallbacks | 2 | 1.17x | 1.14x |
| Mixed row-local fallbacks | 4 | 1.18x | 1.14x |
| Mixed row-local fallbacks | 8 | 1.17x | 1.14x |
| GT-first reordered negative | 0 | 1.21x | 1.22x |
| GT-first reordered negative | 2 | 1.25x | 1.19x |
| GT-first reordered negative | 4 | 1.26x | 1.19x |
| GT-first reordered negative | 8 | 1.22x | 1.18x |
| Two-string float negative | 0 | 0.96x | 0.98x |
| Two-string float negative | 2 | 1.00x | 0.99x |
| Two-string float negative | 4 | 0.99x | 0.98x |
| Two-string float negative | 8 | 1.03x | 1.01x |

## Interpretation

The dynamic path gives a large production-visible win for sample-rich GT-only
VCFs.  On likelihood-heavy rows, it is consistently faster but still limited by
generic per-op work, string/width handling, and IO/compression costs.  Some
float/string-heavy layouts remain near parity or slightly slower than baseline.

## Remaining Work

- Reduce per-sample opcode dispatch in hot FORMAT layouts.
- Improve string and measured-width handling without losing byte identity.
- Consider a later executor-generation layer if generic per-op dispatch remains
  the main gap to historical exact-kernel speed.
