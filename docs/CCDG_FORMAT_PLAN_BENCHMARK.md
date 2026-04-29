# CCDG FORMAT Plan Benchmark Checkpoint

Date: 2026-04-29

Branch: `codex/vcf-avx-sanity`

This file is a checkpoint for the CCDG-oriented FORMAT parser work.  Earlier
versions of the branch used handwritten exact CCDG kernels; those kernels have
now been removed.  The current production candidate is dynamic-only.

## Current Takeaway

The dynamic FORMAT planner is byte-correct on the CCDG subset and larger FORMAT
benchmark corpus.  It is no longer a narrow full-string kernel: it compiles each
FORMAT tag from header metadata and uses one composable executor for supported
tags.

On the CCDG 10k subset, the dynamic-only path is faster than baseline but slower
than the historical handwritten exact kernels.  On GT-only and several
reordered/synthetic workloads, the dynamic path is much closer to the previous
target and can be materially faster than baseline.

## Modes

```sh
HTS_VCF_FORMAT_PLAN=0       # production parser
HTS_VCF_FORMAT_PLAN=1       # dynamic per-tag planner, then production fallback
HTS_VCF_FORMAT_PLAN=interp  # same dynamic planner; manual alias
HTS_VCF_FORMAT_PLAN_STATS=1 # print planner counters from test/test_view
```

## CCDG Data

Source file:

```text
/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/data/original/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

The 10k subset contains 10,000 variant records and 3,202 samples.  The observed
FORMAT distribution is:

| Records | FORMAT |
|---:|---|
| 4,681 | `GT:AB:AD:DP:GQ:PL` |
| 3,774 | `GT:AB:AD:DP:GQ:PGT:PID:PL` |
| 813 | `GT:AD:DP:GQ:PL` |
| 732 | `GT:AD:DP:GQ:PGT:PID:PL` |

The current dynamic planner can compile these layouts from tag metadata rather
than matching the whole FORMAT string.

## Latest Large-Corpus Result

The most recent post-trim run used:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-dynamic-trim-plan \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All planned outputs compared byte-identical to baseline.

| Input | Baseline user | Dynamic plan user | Hits/fallback |
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

## Historical Note

The removed exact kernels remain useful as a performance reference in old
benchmark logs, but they are no longer live code.  New optimization work should
measure `HTS_VCF_FORMAT_PLAN=1` against `HTS_VCF_FORMAT_PLAN=0`; `interp` is
only an alias for manual debugging.

## bcftools Production-Style Check

A clean bcftools `develop` worktree was built against this htslib branch and
timed with:

```sh
bcftools view --no-version -Ob -l 0 [--threads N]
```

on the representative threaded manifest.  All planned outputs compared
byte-identical to baseline.

| Input | Threads | Baseline real | Plan real | Speedup |
|---|---:|---:|---:|---:|
| 1000G chr22 full GT | 0 | 27.48 s | 8.99 s | 3.06x |
| 1000G chr22 full GT | 2 | 26.59 s | 6.99 s | 3.80x |
| 1000G chr22 full GT | 4 | 26.71 s | 6.94 s | 3.85x |
| 1000G chr22 full GT | 8 | 26.62 s | 6.96 s | 3.82x |
| Large CCDG-like synthetic | 0 | 4.43 s | 3.94 s | 1.12x |
| Large CCDG-like synthetic | 2 | 3.46 s | 3.01 s | 1.15x |
| Large CCDG-like synthetic | 4 | 3.47 s | 3.02 s | 1.15x |
| Large CCDG-like synthetic | 8 | 3.46 s | 3.00 s | 1.15x |

## Next Work

- Reduce the CCDG fallback rate without introducing full-string special cases.
- Add selected-sample support so `keep_samples` does not force production
  fallback.
- Lower per-op dispatch and scratch-buffer overhead on likelihood-shaped rows.
- Keep expanding edge fixtures when a new supported FORMAT tag or width pattern
  is added.
