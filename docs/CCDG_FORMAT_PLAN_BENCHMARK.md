# CCDG FORMAT Plan MVP Benchmark

Date: 2026-04-28

Worktree: `/tmp/htslib-vcf-avx-sanity`

Branch: `codex/vcf-avx-sanity`

## Goal

Estimate whether a runtime-planned VCF FORMAT parser can improve end-to-end
compressed VCF/BCF conversion performance on a wide CCDG VCF.

The MVP implementation is gated by:

```sh
HTS_VCF_FORMAT_PLAN=1
HTS_VCF_FORMAT_PLAN_STATS=1
```

It dynamically caches observed FORMAT layouts. The current MVP has direct
executors for the four dominant CCDG layouts:

```text
GT:AB:AD:DP:GQ:PL
GT:AD:DP:GQ:PL
GT:AB:AD:DP:GQ:PGT:PID:PL
GT:AD:DP:GQ:PGT:PID:PL
```

Other layouts fall back to the existing generic FORMAT parser.

## Data

Source file:

```text
/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/data/original/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

Subset used for this benchmark:

```text
/tmp/ccdg_chr22_10k.vcf
```

The subset contains 10,000 variant records plus header lines. It is wide:
3,202 samples and about 866 MiB uncompressed.

Compressed inputs prepared from the subset:

```text
/tmp/ccdg_chr22_10k.vcf.gz
/tmp/ccdg_chr22_10k.bcf
```

Approximate input sizes:

```text
ccdg_chr22_10k.vcf.gz  118 MiB by ls, 129 MiB by du
ccdg_chr22_10k.bcf     152 MiB by ls, 160 MiB by du
```

## FORMAT Coverage

On the 10k CCDG subset after adding `PGT:PID` support:

```text
attempts=10000
hits=10000
fallback=0
parsed_samples=32020000
```

The planned parser therefore handled 100% of records and parsed 32.0 million
sample FORMAT entries directly.

For comparison, before `PGT:PID` support, coverage was:

```text
attempts=10000
hits=5494
fallback=4506
parsed_samples=17591788
```

The fallback records were almost entirely the two layouts containing
`PGT:PID`.

## Four-Cell Compressed Conversion Benchmark

All cells are compressed input to compressed output. Each timing is a single
wall-clock run using `/usr/bin/time -p`; treat these as directional, not a
statistically rigorous benchmark.

| Conversion | Baseline real | FORMAT plan real | Change |
|---|---:|---:|---:|
| VCF.gz -> BCF | 9.150 s | 8.266 s | 9.7% faster |
| BCF -> BCF | 7.168 s | 7.221 s | neutral, 0.7% slower |
| BCF -> VCF.gz | 11.367 s | 11.487 s | neutral, 1.1% slower |
| VCF.gz -> VCF.gz | 13.405 s | 12.670 s | 5.5% faster |

Command shapes:

```sh
./test/test_view.baseline -b -p /tmp/bench_base_vcf_to_bcf.bcf /tmp/ccdg_chr22_10k.vcf.gz
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 ./test/test_view -b -p /tmp/bench_plan_vcf_to_bcf.bcf /tmp/ccdg_chr22_10k.vcf.gz

./test/test_view.baseline -b -p /tmp/bench_base_bcf_to_bcf.bcf /tmp/ccdg_chr22_10k.bcf
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 ./test/test_view -b -p /tmp/bench_plan_bcf_to_bcf.bcf /tmp/ccdg_chr22_10k.bcf

./test/test_view.baseline -z -p /tmp/bench_base_bcf_to_vcf.vcf.gz /tmp/ccdg_chr22_10k.bcf
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 ./test/test_view -z -p /tmp/bench_plan_bcf_to_vcf.vcf.gz /tmp/ccdg_chr22_10k.bcf

./test/test_view.baseline -z -p /tmp/bench_base_vcf_to_vcf.vcf.gz /tmp/ccdg_chr22_10k.vcf.gz
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 ./test/test_view -z -p /tmp/bench_plan_vcf_to_vcf.vcf.gz /tmp/ccdg_chr22_10k.vcf.gz
```

For each cell, the baseline output and planned-parser output were compared with
`cmp` and matched byte-for-byte. The BCF-input cells have
`attempts=0 hits=0 fallback=0` because they never enter the VCF text FORMAT
parser.

## Compressed VCF to Uncompressed BCF

This additional case keeps compressed VCF input but removes output compression
by writing BCF at compression level 0.

| Conversion | Baseline real | FORMAT plan real | Change |
|---|---:|---:|---:|
| VCF.gz -> uncompressed BCF | 2.817 s | 1.930 s | 31.5% faster |

Command shape:

```sh
./test/test_view.baseline -b -l 0 -p /tmp/bench_base_vcfgz_to_ubcf.bcf /tmp/ccdg_chr22_10k.vcf.gz
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 ./test/test_view -b -l 0 -p /tmp/bench_plan_vcfgz_to_ubcf.bcf /tmp/ccdg_chr22_10k.vcf.gz
```

The baseline and planned-parser outputs were compared with `cmp` and matched
byte-for-byte.

## Parse-Only Reference Timings

For context, earlier parse-only tests on the same subsets showed a much larger
effect because output compression was removed from the critical path:

| Dataset | Baseline parse-only | FORMAT plan parse-only | Change |
|---|---:|---:|---:|
| 10k CCDG subset, pre-`PGT:PID` executor | about 2.30 s | about 1.64 s | about 29% faster |
| 100k CCDG subset, pre-`PGT:PID` executor | 23.94 s | 17.71 s | about 26% faster |
| 100k CCDG subset, all-hit executor | 24.22 s | 14.95 s | about 38% faster |
| 100k CCDG VCF.gz -> uncompressed BCF, all-hit executor | 26.65 s | 18.12 s | about 32% faster |

The all-hit executor was byte-identical against baseline on the 10k BCF output
and on a targeted one-record phased-layout test.

## Profiling Notes

After `PGT:PID` support, the generic FORMAT fallback is no longer a meaningful
cost for the CCDG benchmark. A macOS `sample` profile of
`VCF.gz -> uncompressed BCF` on the 100k subset showed the next hot areas:

```text
vcf_plan_parse_int_vector     189 samples
libdeflate input decompress   158 samples
vcf_parse_format              154 samples
bcf_enc_vint                   83 samples
vcf_plan_int_value             42 samples
vcf_plan_copy_string           33 samples
vcf_plan_gt2                   27 samples
vcf_plan_float_value           24 samples
read                           16 samples
```

This is a statistical sample, not exact cycle accounting, but it is useful
directionally. The next parser-side targets are direct integer-vector parsing
for AD/PL and reducing repeated `bcf_enc_vint` work in the planned path.

## Findings

The planned FORMAT parser is viable. With all four dominant CCDG layouts covered,
parse-heavy VCF to uncompressed BCF conversion improves by about 30-40% on the
100k subset.

For fully compressed-to-compressed conversion, output/input compression and VCF
formatting absorb much of the parser win. The MVP still improved VCF-input
conversions by about 5-10%, while BCF-input conversions were unchanged as
expected. When output compression is removed, VCF.gz to uncompressed BCF improves
by about 32%, much closer to the parse-only gain.

The practical takeaway is that FORMAT planning is a better optimization target
than top-level VCF delimiter SIMD scanning. The earlier delimiter-only probe had
100% record coverage but was essentially neutral, while FORMAT planning moved
the parse-heavy workload substantially.

The next highest-value extension is not more FORMAT layout coverage for this
CCDG benchmark, because coverage is already 100%. It is reducing the cost inside
the planned path: AD/PL integer-vector parsing, BCF integer encoding, and then
possibly pipelining decompression/parse/encode once the single-threaded parser
work has been squeezed further.
