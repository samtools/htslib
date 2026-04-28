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

## Follow-Up: Fixed-Width AD/PL Parsing

The first follow-up optimization added fixed-width planned parsers for the most
common biallelic case:

```text
AD width = 2
PL width = 3
```

On the 10k subset, about 82% of records are biallelic, so this removes a large
number of generic integer-vector loop iterations and helper calls while leaving
multi-allelic rows on the generic planned-vector parser.

Correctness checks remained byte-identical against baseline for:

```text
/tmp/ccdg_one_phase.vcf -> uncompressed BCF
/tmp/ccdg_chr22_10k.vcf -> uncompressed BCF
```

Directional timings after the fixed-width parser change:

| Dataset | Previous all-hit plan | Fixed-width AD/PL plan | Change |
|---|---:|---:|---:|
| 100k CCDG VCF -> uncompressed BCF | 14.95 s | 13.1-13.6 s | about 9-12% faster |
| 100k CCDG VCF.gz -> uncompressed BCF | 18.12 s | 15.6-16.5 s | about 9-14% faster |

An attempted range-tracked replacement for `bcf_enc_vint` was also tested. It
preserved byte identity, but it slowed these same parse-heavy cases, so it was
not kept. The likely issue is that tracking ranges during parse adds enough
per-value work to outweigh skipping `bcf_enc_vint`'s later range scan.

## Follow-Up: Compiled Op Interpreter

A second planned-parser tier was added to test whether a more general compiled
FORMAT op interpreter can recover much of the exact-kernel benefit while
covering more layouts.  The exact CCDG kernels still run for
`HTS_VCF_FORMAT_PLAN=1`; `HTS_VCF_FORMAT_PLAN=interp` skips those kernels and
uses only the compiled op interpreter.

Correctness checks:

```text
./test/test_format_plan.sh
HTS_VCF_FORMAT_PLAN=interp ./test/test_view -b -l 0 test/format-plan-edge.vcf
/tmp/ccdg_chr22_10k.vcf -> uncompressed BCF
```

All planned outputs were compared against baseline with `cmp` and matched
byte-for-byte.  The CCDG 10k VCF-input cases had 10,000 attempts, 10,000 hits,
0 fallback, and 32,020,000 parsed samples for both exact and interpreter modes.

Single-pass 10k CCDG conversion matrix, real seconds:

| Conversion | Baseline | Exact kernels | Compiled interp | Exact vs baseline | Interp vs baseline |
|---|---:|---:|---:|---:|---:|
| VCF.gz -> BCF.gz | 9.11 | 7.97 | 8.79 | 12.5% faster | 3.5% faster |
| BCF -> BCF.gz | 7.03 | 7.06 | 7.03 | neutral | neutral |
| BCF -> VCF.gz | 11.20 | 11.32 | 11.21 | neutral | neutral |
| VCF.gz -> VCF.gz | 13.18 | 12.01 | 12.92 | 8.9% faster | 2.0% faster |
| VCF.gz -> uncompressed BCF | 2.79 | 1.64 | 2.61 | 41.2% faster | 6.5% faster |

Parse-heavy uncompressed reference:

| Conversion | Baseline | Exact kernels | Compiled interp | Exact vs baseline | Interp vs baseline |
|---|---:|---:|---:|---:|---:|
| VCF -> uncompressed BCF | 2.56 | 1.36 | 2.33 | 46.9% faster | 9.0% faster |

The compiled interpreter is useful for validating the architecture, but it is
not yet where the performance is.  Its per-sample dynamic dispatch, generic
width pass, generic vector loops, and indirect per-op buffer handling leave it
much closer to the baseline parser than to the exact CCDG kernels.  This argues
for a hybrid approach: use the interpreter as a safe coverage layer and add
small specialized op handlers for the very common shapes inside it, especially
diploid GT, scalar ints, biallelic AD, biallelic PL, and fixed-width strings.

## Follow-Up: Opcode Tape Specialization

The compiled interpreter was then changed from "inspect each op type while
parsing" to a row-specific opcode tape.  The FORMAT string is still cached as a
flexible op list, but after the row width pass each op is resolved to a narrower
handler:

```text
GT2, GT-dynamic, INT1, INT2, INT3, INTN, FLOAT1, FLOATN, STR
```

This preserves the flexible interpreter path for arbitrary defined
String/Integer/Float FORMAT layouts, while avoiding repeated `is_gt` / type
checks and using the same fixed-width integer helpers as the exact CCDG kernel
when the observed row width permits it.

Correctness checks remained byte-identical against baseline for:

```text
./test/test_format_plan.sh
HTS_VCF_FORMAT_PLAN=interp ./test/test_view -b -l 0 test/format-plan-edge.vcf
/tmp/ccdg_chr22_10k.vcf -> uncompressed BCF
```

Single-pass 10k CCDG conversion matrix after opcode specialization, real
seconds:

| Conversion | Baseline | Exact kernels | Opcode interp | Exact vs baseline | Interp vs baseline |
|---|---:|---:|---:|---:|---:|
| VCF.gz -> BCF.gz | 9.19 | 7.99 | 9.28 | 13.1% faster | neutral/noisy |
| BCF -> BCF.gz | 8.04 | 8.22 | 8.10 | neutral | neutral |
| BCF -> VCF.gz | 12.71 | 12.04 | 12.99 | neutral/noisy | neutral/noisy |
| VCF.gz -> VCF.gz | 13.76 | 12.33 | 13.88 | 10.4% faster | neutral/noisy |
| VCF.gz -> uncompressed BCF | 2.87 | 1.68 | 2.43 | 41.5% faster | 15.3% faster |

Parse-heavy uncompressed reference:

| Conversion | Baseline | Exact kernels | Opcode interp | Exact vs baseline | Interp vs baseline |
|---|---:|---:|---:|---:|---:|
| VCF -> uncompressed BCF | 2.57 | 1.42 | 2.12 | 44.7% faster | 17.5% faster |

Relative to the first compiled interpreter measurement, opcode specialization
improved the parse-heavy uncompressed case from 2.33 s to 2.12 s and VCF.gz to
uncompressed BCF from 2.61 s to 2.43 s.  That is real movement, but the exact
kernels remain substantially faster because they also avoid the generic width
measurement, per-op buffer indirection, and per-sample opcode switch.

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
the planned path and then possibly pipelining decompression/parse/encode once
the single-threaded parser work has been squeezed further. After the fixed-width
AD/PL parser, `bcf_enc_vint` and input decompression remain the most obvious
next bottlenecks.
