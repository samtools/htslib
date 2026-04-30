# VCF FORMAT Shape Benchmark Corpus

This directory is a local test and benchmark corpus for the experimental VCF
FORMAT planner in `vcf.c`.  It is intentionally kept under the repository
worktree instead of `/tmp` so the inputs survive restarts.

The canonical feature docs are:

- `docs/FORMAT_PLAN_OVERVIEW.md` for the high-level feature summary;
- `docs/FORMAT_PLAN_CURRENT.md` for the current implementation and benchmark
  tables;
- `docs/FORMAT_PLAN_EXPERIMENT_LOG.md` for the historical experiment log.

## Layout

```text
bench/format-shape/
  inputs.tsv                 input manifest used by the benchmark script
  public/                    downloaded public VCF slices
  synthetic/                 generated VCFs covering targeted FORMAT shapes
  large/                     meaningful multi-second benchmark inputs/results
  scripts/make_synthetic.pl  deterministic synthetic VCF generator
  scripts/make_large_synthetic.pl
  scripts/run_bench.sh       baseline/plan timing and cmp runner
  scripts/run_thread_bench.sh threaded timing and cmp runner
  scripts/run_bcftools_bench.sh bcftools threaded timing runner
  scripts/run_bcftools_command_bench.sh broader bcftools command runner
  scripts/run_bcftools_command_bench_stream.sh checksum-only large-output runner
  results/                   generated timing logs and BCF outputs
```

The downloaded/generated VCF inputs and benchmark result files are intentionally
ignored by git to avoid accidentally pushing large local data.  The manifests,
scripts, and docs are tracked; local data and timing summaries can be
regenerated from the commands below.

`results/` can be regenerated at any time and may become large.  The scripts
write timing/check summaries and keep BCF outputs locally so `cmp` checks are
inspectable, but `.gitignore` excludes those rerun artifacts.

## Repo Tests

The small correctness cases that should travel with the implementation now live
in the normal htslib test harness, not only in this benchmark directory.
`make check` runs `test_vcf_format_plan` inside `test/test.pl` plus
`test/test_format_plan_cache`.  Those tests assert byte-identical planned output
at the parser-output level, selected-sample behavior, rollback after partial
planned parsing, malformed-input failure behavior, and header-cache generation
invalidation.  Fallback reason counters remain local diagnostics for benchmark
analysis rather than production test assertions.

The benchmark corpus remains for performance and production-shape coverage.  It
should not become a normal test-suite dependency because several inputs are
large public VCFs or generated multi-second workloads.

## Public Inputs

The small `public/` and `synthetic/` inputs are smoke/correctness fixtures.  They
are not large enough to provide stable timing signal except for the CCDG 10k
subset.  Use `large/inputs.tsv` for optimization decisions.

The public files were sliced with `tabix -h URL REGION | ./bgzip -c > file`.
They are small enough to keep in the worktree but diverse enough to catch
non-FORMAT and real-world INFO-heavy workloads.

| File | Source | Shape |
|---|---|---|
| `public/ccdg_chr22_10k.vcf.gz` | local CCDG subset | 3,202-sample CCDG likelihood FORMAT |
| `public/1000g_chr22_genotypes_16050k_16150k.vcf.gz` | 1000 Genomes Phase 3 chr22 genotypes | sample-rich `GT` FORMAT |
| `public/1000g_wgs_sites_chr22_16050k_16300k.vcf.gz` | 1000 Genomes Phase 3 WGS sites | sites-only |
| `public/clinvar_grch38_chr22_16050k_20000k.vcf.gz` | ClinVar GRCh38 VCF | sites-only clinical annotations |
| `public/gnomad_v4.1_exomes_sites_chr22_20000k_20100k.vcf.gz` | gnomAD v4.1 exomes chr22 | sites-only, INFO-heavy |
| `large/public/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` | GIAB HG002 v4.2.1 | 4,048,342-record single-sample truth-set small variants |
| `large/public/giab/HG002_GRCh38_v5.0q_smvar.vcf.gz` | GIAB HG002 v5.0q GRCh38 | 5,945,525-record single-sample small variants |
| `large/public/giab/HG002_GRCh38_v5.0q_stvar.vcf.gz` | GIAB HG002 v5.0q GRCh38 | 6,268,852-record single-sample structural variants |
| `large/public/giab/HG002_CHM13v2.0_v5.0q_smvar.vcf.gz` | GIAB HG002 v5.0q CHM13v2.0 | 5,829,374-record single-sample small variants |

Source URLs used:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_smvar.vcf.gz
```

The parent CCDG/1000G high-coverage chr22 file for
`public/ccdg_chr22_10k.vcf.gz` is:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

It is 26.0 GiB compressed.  For local reruns, point the full-CCDG manifest at a
local copy such as:

```text
/path/to/local/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

Do not run the normal output-materializing command harness on this file.  A
single uncompressed BCF output reached 155 GiB before the run was interrupted.
Use the streaming checksum harness below instead.

## Synthetic Inputs

The synthetic files are generated by:

```sh
bench/format-shape/scripts/make_synthetic.pl bench/format-shape/synthetic
for f in bench/format-shape/synthetic/*.vcf; do ./bgzip -f "$f"; done
```

They cover:

- CCDG-like likelihood layouts with optional `AB` and `PGT/PID`,
- reordered likelihood fields,
- fixed numeric vectors,
- float-vector plus string FORMAT fields,
- multiallelic AD/PL likelihood rows.

## Running

Build the tools first:

```sh
make test/test_view tabix bgzip
```

Run all inputs:

```sh
bench/format-shape/scripts/run_bench.sh
```

Run only the meaningful large corpus:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results \
  bench/format-shape/scripts/run_bench.sh bench/format-shape/large/inputs.tsv
```

`KEEP_OUTPUTS=0` still writes temporary BCF files and compares them with `cmp`,
but deletes the large BCF outputs after each input is checked.

Run the threaded scaling corpus:

```sh
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-threaded \
  bench/format-shape/scripts/run_thread_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

By default this runs unthreaded plus `-@ 2`, `-@ 4`, and `-@ 8`.  Override with
`THREADS_LIST="2 4 8"` or a similar space-separated list.  The threaded manifest
now mirrors the full large corpus so thread scaling is checked across the same
real and synthetic workload shapes as the primary benchmark.

The script runs each input in two modes:

```text
baseline: HTS_VCF_FORMAT_PLAN=0
plan:     HTS_VCF_FORMAT_PLAN=1
```

It writes:

```text
bench/format-shape/results/timings.tsv
bench/format-shape/results/checks.tsv
```

`checks.tsv` compares plan BCF output against baseline with `cmp`.
The threaded runner writes the same files under its selected output directory,
with an additional `threads` column.

Run the same threaded corpus through bcftools:

```sh
BCFTOOLS=/path/to/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools \
  bench/format-shape/scripts/run_bcftools_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

This uses `bcftools view --no-version -Ob -l 0`, compares planned output against
baseline with `cmp`, and records the same `0 2 4 8` thread counts by default.
It does not report planner counters because bcftools does not expose the
`test/test_view` stats hook.

To exercise selected-sample parsing, set `SAMPLE_COUNT=N`.  The runner queries
the first N samples from each input with `bcftools query -l` and passes them to
`bcftools view -s`; sites-only inputs have no sample list and run unchanged.

```sh
BCFTOOLS=/path/to/bcftools SAMPLE_COUNT=2 \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-keep2 \
  bench/format-shape/scripts/run_bcftools_bench.sh \
  bench/format-shape/large/threaded-inputs.tsv
```

Run broader bcftools command shapes:

```sh
BCFTOOLS=/path/to/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-commands \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-command-inputs.tsv
```

This runner is intended to be a bridge toward future tests.  It runs each
command once with `HTS_VCF_FORMAT_PLAN=0` and once with
`HTS_VCF_FORMAT_PLAN=1`, then compares outputs with `cmp`.

The default command set is:

| Command | Purpose | Output check |
|---|---|---|
| `view_bcf` | Full `bcftools view --no-version -Ob -l 0` conversion. | Binary BCF `cmp`. |
| `view_sites` | `bcftools view -G` after dropping genotypes. | Binary BCF `cmp`. |
| `query_sites` | Fixed-column query that should not benefit from FORMAT parsing. | Text `cmp`. |
| `query_format` | Query `%GT` for the first `QUERY_SAMPLE_COUNT` samples. | Text `cmp`. |
| `stats` | `bcftools stats` over the input. | Text `cmp`. |
| `filter_gt` | `bcftools view -i 'GT="alt"'` for the first `QUERY_SAMPLE_COUNT` samples. | Binary BCF `cmp`. |
| `merge_self` | `bcftools merge --no-index --force-samples` of the input with itself. | Binary BCF `cmp`. |

`query_format`, `filter_gt`, and `merge_self` are skipped for sites-only inputs.
By default the query/filter commands select two samples
(`QUERY_SAMPLE_COUNT=2`) to avoid generating enormous text output on cohort-scale
VCFs.  Override with:

```sh
COMMANDS="query_format stats" QUERY_SAMPLE_COUNT=8 THREADS_LIST="0 4" \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-command-inputs.tsv
```

The runner writes:

```text
timings.tsv   name, command, threads, mode, real/user/sys
checks.tsv    baseline-vs-plan cmp status, including skipped_no_samples
commands.tsv  command descriptions captured with the result directory
```

For very large inputs, use the streaming checksum variant.  It runs the same
command families but pipes output through `cksum` and compares checksums instead
of storing complete BCF/text outputs:

```sh
BCFTOOLS=/path/to/bcftools \
OUTDIR=bench/format-shape/large/results-bcftools-full-ccdg-stream \
  bash bench/format-shape/scripts/run_bcftools_command_bench_stream.sh \
  bench/format-shape/large/bcftools-full-ccdg-inputs.tsv
```

The full CCDG chr22 streaming run wrote:

```text
bench/format-shape/large/results-bcftools-full-ccdg-stream/timings.tsv
bench/format-shape/large/results-bcftools-full-ccdg-stream/checks.tsv
bench/format-shape/large/results-bcftools-full-ccdg-stream/checksums.tsv
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

For CI, the likely future shape is to keep one or two tiny inputs per command
and assert `checks.tsv` has only `ok` or expected `skipped_no_samples` rows.
The large corpus should remain a performance benchmark rather than a normal
test-suite dependency.

Run the GIAB plus CCDG correctness/performance pass:

```sh
BCFTOOLS=/path/to/bcftools \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-giab-ccdg-prod-hardening \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-giab-ccdg-inputs.tsv
```

If using the sibling bcftools checkout in this workspace, build it against this
HTSlib checkout explicitly:

```sh
cd ../bcftools-htslib-vcf-plan
make HTSDIR=../htslib-vcf-avx-sanity bcftools
```

This pass is primarily a production-shape correctness check.  GIAB is
single-sample, so it does not show the large cohort speedups, but it does cover
real truth-set small-variant and structural-variant FORMAT details.  The first
GIAB v5.0q run exposed a planned-path bug where `.|.` was serialized as `./.`.
The GT2 parser now preserves phased missing alleles, and the fixed rerun has
all baseline-vs-plan command outputs comparing `ok`.

Latest hardened GIAB/CCDG command run:

| Input | Command | Real speedup | User speedup |
|---|---|---:|---:|
| CCDG 10k | view_bcf | 1.14x | 1.14x |
| CCDG 10k | view_sites | 1.13x | 1.14x |
| CCDG 10k | query_format | 1.52x | 1.56x |
| CCDG 10k | filter_gt | 1.12x | 1.12x |
| GIAB HG002 GRCh38 v4.2.1 | view_bcf | 1.09x | 1.09x |
| GIAB HG002 GRCh38 v4.2.1 | query_format | 1.07x | 1.07x |
| GIAB HG002 GRCh38 v4.2.1 | filter_gt | 1.09x | 1.09x |
| GIAB HG002 GRCh38 v5.0q small variants | view_bcf | 1.09x | 1.09x |
| GIAB HG002 GRCh38 v5.0q small variants | query_format | 1.09x | 1.07x |
| GIAB HG002 GRCh38 v5.0q structural variants | view_bcf | 1.09x | 1.09x |
| GIAB HG002 GRCh38 v5.0q structural variants | query_format | 1.02x | 1.02x |
| GIAB HG002 CHM13 v5.0q small variants | view_bcf | 1.07x | 1.07x |
| GIAB HG002 CHM13 v5.0q small variants | query_format | 1.06x | 1.06x |

`merge_self` is intentionally not in the default `COMMANDS` list because it can
produce very large outputs on cohort-scale inputs.  Run it against the smaller
merge manifest:

```sh
BCFTOOLS=/path/to/bcftools COMMANDS=merge_self \
KEEP_OUTPUTS=0 OUTDIR=bench/format-shape/large/results-bcftools-merge \
  bench/format-shape/scripts/run_bcftools_command_bench.sh \
  bench/format-shape/large/bcftools-merge-inputs.tsv
```

This is not a semantic recommendation to merge a file with itself in production;
it is a controlled benchmark shape.  `--force-samples` creates distinct sample
names and `--no-index` avoids needing local tabix indexes for generated slices.

The latest local merge run wrote:

```text
bench/format-shape/large/results-bcftools-merge/timings.tsv
bench/format-shape/large/results-bcftools-merge/checks.tsv
```

All planned merge outputs compared byte-identical to baseline.  The small 1000G
genotype input improved from 0.14 s to 0.10 s, the 1024-sample CCDG-like input
improved from 4.50 s to 4.33 s, and the 1024-sample float/string input was
unchanged at 2.69 s.

## Large Corpus

`large/inputs.tsv` currently contains:

- the CCDG 10k subset,
- the full 1000 Genomes chr22 Phase 3 genotype VCF,
- eight generated 2,048-sample synthetic FORMAT workloads:
  CCDG-like likelihood, reordered likelihood, multiallelic likelihood,
  float/string FORMAT, variable phase-string widths, row-local likelihood
  fallbacks, GT-first wrong-order likelihood-like rows, and two-string
  float rows.

`large/threaded-inputs.tsv` mirrors this full corpus for `-@` scaling checks.
`large/bcftools-command-inputs.tsv` is a smaller representative set for the
broader command benchmark: GT-only, real CCDG-like FORMAT, reordered FORMAT,
string/float negative control, and an INFO-heavy sites-only gnomAD slice.
`large/bcftools-merge-inputs.tsv` is smaller still, so merge output does not
explode during routine local benchmarks.

To refresh only the newer cache-regression synthetic files without rewriting the
older large VCFs:

```sh
SYNTHETIC_ONLY_NEW=1 \
  bench/format-shape/scripts/make_large_synthetic.pl \
  bench/format-shape/large/synthetic 2048
```

The latest large run used this local output directory:

```text
bench/format-shape/large/results-prod-hardening2/timings.tsv
bench/format-shape/large/results-prod-hardening2/checks.tsv
```

Generated result files are ignored; the summary below is the portable record.
All plan outputs in that run compared byte-identical to baseline.

That run includes fallback reason diagnostics.  In the CCDG 10k slice, the
planner hit 9,861 of 10,000 rows; the remaining 139 rows fell back for
`string_width`, meaning their measured string field exceeded the current
256-byte planned cap.

One rejected optimization is recorded in
`bench/format-shape/large/results-opt-nosubset-split`: splitting the all-samples
loop from the `keep_samples` loop preserved correctness but slowed the planned
rows, so that code was reverted.
