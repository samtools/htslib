# Dynamic FORMAT Plan: Current Implementation

This document describes the implementation currently present in `vcf.c`, the
correctness boundaries, and the latest benchmark results.

## Entry Point

`vcf_parse_format()` first calls `vcf_parse_format_planned()` when
`HTS_VCF_FORMAT_PLAN` is enabled.  The planned path either parses the whole
FORMAT column or returns `-3`, allowing the existing generic parser to handle the
column unchanged.

```text
HTS_VCF_FORMAT_PLAN enabled
  -> fetch or compile header-owned FORMAT/header plan
  -> resolve row-local widths
  -> composable executor
  -> generic fallback on unsupported or suspicious rows
```

The only enabled spelling is `HTS_VCF_FORMAT_PLAN=1`.  Unknown values are
treated as disabled so typos such as `off` or `false` do not accidentally enable
the planner.

## Plan Compilation

Plans are cached in private `bcf_hdr_aux_t` state by literal FORMAT string plus
the header's private FORMAT-plan generation.  This is important because VCF
header IDs, declared types, and number models are header-local.  The cache grows
from 16 entries up to 128 entries, uses heap storage for long FORMAT strings,
and also caches unsupported schemas so repeated odd rows do not repeatedly pay
compile cost.

Planner statistics are collected only when `HTS_VCF_FORMAT_PLAN_STATS=1` is
also set.  Normal production parsing therefore avoids touching the process-wide
test counters.  The test hook reports both aggregate attempts/hits/fallbacks and
fallback reason counters: unsupported schema, numeric width, string width, GT
shape, parse failure, separator mismatch, and sample-count mismatch.

`bcf_hdr_sync()` clears the header-owned plan cache and increments the private
generation after header dictionaries are rebuilt.  The planner also refuses to
compile while `h->dirty` is set, leaving unsynced or header-repair cases on the
generic parser.

The cache is mutable header-owned state, like other htslib header scratch
storage.  Callers should not concurrently parse through the same `bcf_hdr_t`
from multiple threads.

The compile step rejects:

- undefined FORMAT tags;
- duplicate FORMAT tags;
- unsupported header types;
- unsupported number models;
- `GT` declarations that are not `Type=String,Number=1`.
- measured-string plus float-vector schemas that do not also have integer-vector
  work for the planned executor.

Undefined tags intentionally fall back to the generic parser so existing
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
fields use a row-local measurement pass.  Numeric vectors remain capped at 64
values per FORMAT field in the planned path.  Measured strings are capped
separately at 256 bytes per row field, which keeps common phase-set annotations
on the fast path while bounding scratch-buffer and transposition work.

## Executor

BCF stores FORMAT data transposed by tag: all samples for FORMAT op 0, then all
samples for FORMAT op 1, and so on.  The dynamic executor parses VCF samples in
sample-major order and writes that transposed BCF layout.

Leading fixed-width `GT2` and `FLOAT1` rows can write directly into `v->indiv`.
Other rows are staged in header scratch memory, then encoded after sample
parsing so integer range and observed-width metadata are known.

For fixed-width vector fields, the executor can compact underfilled rows to the
observed row maximum before BCF encoding.  This avoids whole-row fallback when
the generic parser would also emit a narrower byte-identical vector width.

## Fallback Policy

Supported cached plans are probed on every row.  If row-local validation fails,
the executor rolls back its partial `v->indiv` writes and the generic parser
handles the whole FORMAT column for that record.  The fallback does not disable
or cool down the cached plan; nearby rows with the same FORMAT schema can still
take the optimized path.

Compile-time unsupported schemas are still cached as unsupported, so repeated
unoptimizable FORMAT strings pay the compile/classification cost once and then
fall back directly to the generic parser.

## Correctness Rules

The planned parser must preserve these invariants:

- no planned parsing while the header has unsynced dictionary changes;
- header IDs, types, and number models are resolved before execution;
- selected-sample parsing must honor `h->keep_samples`, use `h->nsamples_ori`
  for input-column scans, and set `v->n_sample` to the retained sample count;
- duplicate or undefined tags use the generic parser;
- measured-string plus float-vector schemas without integer-vector work use the
  generic parser;
- unsupported GT encodings force fallback;
- numeric vectors preserve observed width and vector-end padding;
- strings use observed maximum byte length and zero-pad shorter samples;
- integer and float overflow/error behavior must match production htslib or
  force fallback;
- successful planned rows run the same final FORMAT consistency check as the
  generic parser via `vcf_parse_format_check7()`;
- direct writes to `v->indiv` must roll back before fallback.

Focused validation lives in the existing `test/test.pl` harness as
`test_vcf_format_plan`.  It compares generic parsing and
`HTS_VCF_FORMAT_PLAN=1` byte-for-byte with `cmp`, and also verifies that
unrecognized control values such as `HTS_VCF_FORMAT_PLAN=off` behave like the
generic parser.  The repo fixtures cover numeric-width and GT-shape fallback,
mixed string/float schemas kept on the generic parser, cache growth, long FORMAT
strings, string-width fallback, separator fallback, parse fallback with rollback,
repeated wide GT values, float-vector compaction, selected-sample skipping of
malformed unselected columns, and sample-count mismatch.  The selected-sample
checks compare explicit inclusion and exclusion lists (`S1,S3`, `S2`, and
`^S2`) and also verify retained-sample float widths do not depend on skipped
input columns.
`test/test_format_plan_cache` mutates and resyncs a header after a plan has been
compiled for the same FORMAT string, then verifies the row is planned again with
the new metadata.

## Large Corpus Benchmark

Command:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-prod-hardening2 \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

All planned outputs compared byte-identical to baseline.

| Input | Baseline user | Plan user | User speedup | Hits/fallback |
|---|---:|---:|---:|---:|
| CCDG 10k | 2.47 s | 2.21 s | 1.12x | 9,861 / 139 |
| 1000G chr22 full GT | 24.61 s | 9.48 s | 2.60x | 1,103,547 / 0 |
| Large CCDG-like synthetic | 4.00 s | 3.68 s | 1.09x | 20,000 / 0 |
| Large reordered likelihood | 2.86 s | 2.42 s | 1.18x | 20,000 / 0 |
| Large multiallelic likelihood | 3.08 s | 2.67 s | 1.15x | 16,000 / 0 |
| Large float/string | 2.88 s | 2.86 s | 1.01x | 0 / 16,000 |
| Variable phase widths | 2.53 s | 2.45 s | 1.03x | 12,000 / 0 |
| Mixed row-local fallbacks | 2.14 s | 1.84 s | 1.16x | 12,000 / 0 |
| GT-first reordered | 1.68 s | 1.41 s | 1.19x | 12,000 / 0 |
| Two-string float | 2.20 s | 2.19 s | 1.00x | 0 / 12,000 |

The CCDG 10k fallbacks are all `string_width=139`, meaning only rows with
measured string fields wider than the 256-byte planned cap use the generic
parser.  The float/string control fixtures still fall back as unsupported
because the mixed string/float shape boundary keeps those rows on the generic
parser.  Briefly tested runtime guards regressed sparse-fallback CCDG-like
layouts, so the current implementation leaves row-local fallbacks local to the
record.

## Full Threaded Corpus Benchmark

Command:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-threaded-profit-gate \
  bench/format-shape/scripts/run_thread_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

All 40 planned outputs compared byte-identical to baseline.  Generated result
files are ignored; the table below summarizes the recorded real-time speedup
from `bench/format-shape/large/results-threaded-profit-gate`.

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
/path/to/bcftools-htslib-vcf-plan
```

using this htslib worktree:

```sh
make HTSDIR=/path/to/htslib bcftools
```

Timing command:

```sh
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
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
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
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

## bcftools Command Benchmark

The broader command runner exercises bcftools paths that either consume FORMAT
records, discard FORMAT records, or mostly operate on site-level data:

```sh
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-commands \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-command-inputs.tsv
```

All applicable planned outputs compared byte-identical to baseline.  FORMAT
commands were skipped for the sites-only gnomAD input as expected.

| Input | Command | Real speedup | User speedup |
|---|---|---:|---:|
| CCDG 10k | view_bcf | 1.11x | 1.12x |
| CCDG 10k | view_sites | 1.12x | 1.13x |
| CCDG 10k | query_format | 1.51x | 1.56x |
| CCDG 10k | filter_gt | 1.11x | 1.12x |
| 1000G chr22 full GT | view_bcf | 2.79x | 2.94x |
| 1000G chr22 full GT | view_sites | 2.98x | 3.02x |
| 1000G chr22 full GT | query_format | 1.94x | 1.94x |
| 1000G chr22 full GT | filter_gt | 1.57x | 1.58x |
| Large reordered likelihood | view_bcf | 1.21x | 1.22x |
| Large reordered likelihood | view_sites | 1.20x | 1.20x |
| Large reordered likelihood | query_format | 1.39x | 1.42x |
| Large reordered likelihood | filter_gt | 1.14x | 1.14x |
| Large float/string | view_bcf | 1.02x | 1.02x |
| Large float/string | query_format | 1.01x | 1.00x |
| gnomAD sites chr22 | view_bcf | 0.98x | 1.00x |
| gnomAD sites chr22 | query_sites | 1.00x | 1.08x |

`query_sites` and `stats` were generally neutral because they do little or no
FORMAT work.  The small negative rows, such as CCDG `stats` at 0.94x real and
float/string `stats` at 0.93x real, are still within the area to watch for
planner overhead in workloads that do not benefit from FORMAT decoding.

## bcftools Merge Benchmark

`merge_self` is kept out of the default command list because merge output can
grow quickly.  It was run against the smaller merge manifest:

```sh
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
COMMANDS=merge_self KEEP_OUTPUTS=0 \
OUTDIR=bench/format-shape/large/results-bcftools-merge \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-merge-inputs.tsv
```

All planned merge outputs compared byte-identical to baseline.

| Input | Baseline real | Plan real | Real speedup | Baseline user | Plan user |
|---|---:|---:|---:|---:|---:|
| Small 1000G genotypes | 0.14 s | 0.10 s | 1.40x | 0.13 s | 0.08 s |
| Large CCDG likelihood 1024s | 4.50 s | 4.33 s | 1.04x | 4.05 s | 3.91 s |
| Large float/string 1024s | 2.69 s | 2.69 s | 1.00x | 2.40 s | 2.41 s |

## GIAB and CCDG Command Check

GIAB HG002 files were added as real-world single-sample correctness fixtures:
NIST v4.2.1 GRCh38 small variants, v5.0q GRCh38 small variants, v5.0q GRCh38
structural variants, and v5.0q CHM13v2.0 small variants.  The same bcftools
command suite was run against those files plus the all-sample CCDG 10k slice:

```sh
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-giab-ccdg-prod-hardening \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-giab-ccdg-inputs.tsv
```

The initial GIAB v5.0q run found a correctness bug: the planned GT2 parser
encoded phased missing genotypes such as `.|.` as unphased `./.`.  The parser
now preserves the phase bit for missing alleles, and `test/format-plan-edge.vcf`
has an explicit phased-missing GT row.  After the fix, every command output in
this run compared byte-identical/text-identical to baseline.  The table below
shows user-time speedups from the latest hardened rerun.

| Input | Records | Samples | view_bcf | query_format | filter_gt | Notes |
|---|---:|---:|---:|---:|---:|---|
| CCDG 10k | 10,000 | 3,202 | 1.14x | 1.56x | 1.12x | Cohort FORMAT win remains visible. |
| GIAB HG002 GRCh38 v4.2.1 | 4,048,342 | 1 | 1.09x | 1.07x | 1.09x | Single-sample truth-set small variants. |
| GIAB HG002 GRCh38 v5.0q small variants | 5,945,525 | 1 | 1.09x | 1.07x | 1.07x | Includes phased missing GTs. |
| GIAB HG002 GRCh38 v5.0q structural variants | 6,268,852 | 1 | 1.09x | 1.02x | 1.08x | Structural-variant FORMAT coverage. |
| GIAB HG002 CHM13 v5.0q small variants | 5,829,374 | 1 | 1.07x | 1.06x | 1.16x | Alternate reference truth-set coverage. |

The parent CCDG/1000G high-coverage chr22 file is 26.0 GiB compressed:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

For local reruns, point the full-CCDG manifest at a local copy such as:

```text
/path/to/local/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

The normal command harness is unsafe for this input because one full
`view_bcf -Ob -l 0` output reached 155 GiB before the run was interrupted.  The
full-file benchmark therefore used the streaming checksum harness:

```sh
BCFTOOLS=/path/to/bcftools-htslib-vcf-plan/bcftools \
OUTDIR=bench/format-shape/large/results-bcftools-full-ccdg-stream \
  bash bench/format-shape/scripts/run_bcftools_command_bench_stream.sh \
  bench/format-shape/large/bcftools-full-ccdg-inputs.tsv
```

All baseline-vs-plan checksums compared `ok`.

| Command | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 678.46 s | 562.96 s | 1.21x | 476.41 s | 377.47 s | 1.26x |
| `view_sites` | 472.27 s | 403.28 s | 1.17x | 455.70 s | 386.18 s | 1.18x |
| `query_sites` | 71.44 s | 76.78 s | 0.93x | 67.02 s | 72.00 s | 0.93x |
| `query_format` | 124.14 s | 76.88 s | 1.61x | 119.16 s | 72.27 s | 1.65x |
| `stats` | 77.45 s | 77.12 s | 1.00x | 72.86 s | 72.55 s | 1.00x |
| `filter_gt` | 531.20 s | 453.21 s | 1.17x | 512.95 s | 434.35 s | 1.18x |

## Executor Optimization Pass

The latest optimization pass stayed within the generic per-op executor.  It did
not add schema-specific kernels.  Retained changes are:

- skip observed-count tracking for row ops that cannot compact;
- update integer ranges directly on the common positive-integer path;
- fail over-wide measured fields during the measurement pass;
- remove nullable `nread` checks from planner-private vector helpers.

Focused tests passed:

```sh
make test/test_view test/test_format_plan_cache
cd test && REF_PATH=: ./test.pl -F vcf_format_plan
test/test_format_plan_cache
git diff --check
```

The htslib large corpus result is in
`bench/format-shape/large/results-opt-batch1b`.  All planned outputs compared
byte-identical to baseline.

| Input | Plan user | User speedup | Hits/fallback |
|---|---:|---:|---:|
| CCDG 10k | 2.20 s | 1.14x | 8,396 / 1,604 |
| 1000G chr22 full GT | 8.99 s | 2.79x | 1,103,547 / 0 |
| Large CCDG-like synthetic | 3.68 s | 1.09x | 20,000 / 0 |
| Large reordered likelihood | 2.38 s | 1.22x | 20,000 / 0 |
| Large multiallelic likelihood | 2.64 s | 1.21x | 16,000 / 0 |
| Large float/string | 2.88 s | 1.00x | 0 / 16,000 |
| Variable phase widths | 2.44 s | 1.05x | 12,000 / 0 |
| Mixed row-local fallbacks | 1.83 s | 1.20x | 12,000 / 0 |
| GT-first reordered | 1.41 s | 1.23x | 12,000 / 0 |
| Two-string float | 2.24 s | 1.00x | 0 / 12,000 |

The `keep_samples`/all-samples loop split was tested and rejected.  It preserved
correctness, but `bench/format-shape/large/results-opt-nosubset-split` was
slower across the planned rows, so the change was reverted.

For bcftools-level validation, the sibling bcftools checkout must be built
against this checkout explicitly:

```sh
make HTSDIR=../htslib-vcf-avx-sanity bcftools
```

The standard GIAB/CCDG command result is in
`bench/format-shape/large/results-bcftools-giab-ccdg-opt-batch1`; all outputs
compared `ok`.  CCDG 10k user-time speedups were 1.12x for `view_bcf`, 1.55x
for `query_format`, and 1.11x for `filter_gt`.  GIAB single-sample FORMAT query
rows were roughly 1.08-1.12x faster; site-only controls and `stats` remain
neutral/noisy as expected.

## Fallback Diagnostics And String Width Tuning

A later pass added fallback reason counters and split the planned width cap
into numeric and string limits:

- numeric measured vectors remain capped at 64 values;
- measured strings are capped at 256 bytes;
- numeric/string width fallbacks are counted but do not disable the cached plan.

A 512-byte string cap was tested first.  It recovered all CCDG 10k planner
fallbacks, but the bcftools-level signal was mixed.  The retained 256-byte cap
keeps almost all CCDG rows on the planned path while leaving the longest string
rows on the generic parser.

Focused CCDG 10k htslib result at 256 bytes:

| Metric | Value |
|---|---:|
| Baseline user | 2.43 s |
| Plan user | 2.15 s |
| Hits / fallback | 9,861 / 139 |
| Fallback reason | `string_width=139` |

The standard GIAB/CCDG bcftools command result for the retained version is in
`bench/format-shape/large/results-bcftools-giab-ccdg-cap256`; all outputs
compared `ok`.

| Input | `view_bcf` user | `query_format` user | `filter_gt` user |
|---|---:|---:|---:|
| CCDG 10k | 1.13x | 1.56x | 1.10x |
| GIAB HG002 GRCh38 v4.2.1 | 1.08x | 1.08x | 1.04x |
| GIAB HG002 GRCh38 v5.0q small variants | 1.13x | 1.08x | 1.03x |
| GIAB HG002 GRCh38 v5.0q structural variants | 1.11x | 1.15x | 1.04x |
| GIAB HG002 CHM13 v5.0q small variants | 1.08x | 1.07x | 1.03x |

## Repo Test Harness Hardening

The latest hardening pass moved the important correctness checks into the normal
htslib `test/test.pl` harness instead of leaving them only in
`bench/format-shape`.  The `make check` coverage is intentionally black-box and
includes:

- byte-identity checks for all small planned-path fixtures;
- generic parser vs planned parser comparisons;
- disabled-control comparisons for `HTS_VCF_FORMAT_PLAN=off`;
- a rollback row where planned parsing starts and then falls back after a DP
  overflow;
- repeated unsupported wide GT values;
- selected-sample parsing where malformed unselected sample fields must be
  skipped and must not affect emitted widths;
- malformed sample-count input, where both generic and planned modes must fail;
- cache-generation coverage in `test/test_format_plan_cache`.

The planned executor now calls `vcf_parse_format_check7()` on success, so the
planned path shares the generic parser's final FORMAT cardinality check.  The
fallback counters are test-only diagnostics, exposed through renamed
`*_for_test` hooks rather than API-looking `hts_*` names.

## Interpretation

The dynamic path gives a large production-visible win for sample-rich GT-only
VCFs.  On likelihood-heavy rows, it is consistently faster but still limited by
generic per-op work, string/width handling, and IO/compression costs.  Some
float/string-heavy layouts remain near parity or slightly slower than baseline.
The broader bcftools command run supports the same story: commands that expose
FORMAT parsing benefit; commands dominated by site-only logic, stats, merge
bookkeeping, or compression are neutral.

## Remaining Work

- Reduce per-sample opcode dispatch in hot FORMAT layouts.
- Improve string and measured-width handling without losing byte identity.
- Consider a later executor-generation layer if generic per-op dispatch remains
  the main gap to historical exact-kernel speed.
