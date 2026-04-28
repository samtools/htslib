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

## Follow-Up: Strict Width and Shape Executors

The next iteration hardened correctness and tested more aggressive FORMAT
planning:

- planned integer parsing now detects BCF int32 payload overflow and falls back,
  avoiding undefined overflow and preserving generic warning/missing behavior;
- exact AD/PL paths validate that the observed max vector width matches the
  emitted width, falling back for sparse rows that generic htslib would encode
  narrower;
- the general interpreter can skip the observed-width pass for strict
  header/allele-count-derived numeric rows;
- common numeric opcode tapes `GT2:INT2:INT1:INT1:INT3` and
  `GT2:FLOAT1:INT2:INT1:INT1:INT3` use shape-level executors;
- validated `GT2` rows emit BCF int8 directly instead of calling
  `bcf_enc_vint()`.

`test/format-plan-edge.vcf` now includes an all-missing AD/PL row to verify that
the exact path falls back when its expected vector width would not match generic
observed-width BCF output.

Correctness checks remained byte-identical:

```text
./test/test_format_plan.sh
/tmp/ccdg_chr22_10k.vcf -> uncompressed BCF, exact mode
/tmp/ccdg_chr22_10k.vcf -> uncompressed BCF, interpreter mode
```

Parse-heavy 10k CCDG reference after these changes:

| Conversion | Baseline | Exact kernels | Strict/shape interp |
|---|---:|---:|---:|
| VCF -> uncompressed BCF | 2.63 s | 1.61 s | 2.31 s |

Full 10k CCDG compressed matrix, real seconds:

| Conversion | Baseline | Exact kernels | Strict/shape interp |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz | 9.26 | 8.22 | 8.94 |
| BCF -> BCF.gz | 7.18 | 7.20 | 7.16 |
| BCF -> VCF.gz | 11.45 | 11.33 | 11.85 |
| VCF.gz -> VCF.gz | 14.37 | 13.55 | 13.52 |
| VCF.gz -> uncompressed BCF | 2.94 | 1.92 | 2.66 |

On a 3k CCDG subset containing only non-phase FORMAT layouts, the strict/shape
interpreter improved over baseline but still did not approach the exact kernel:

| Dataset | Baseline | Exact kernels | Strict/shape interp |
|---|---:|---:|---:|
| 3k non-phase VCF -> uncompressed BCF | 0.68 s | 0.34 s | 0.58 s |

The takeaway is mixed.  The hardening is worth keeping, and direct `GT2`
encoding is simple and safe.  However, shape-level dispatch alone does not close
the remaining gap.  The next high-ROI parser-side experiment should reduce
memory traffic by parsing validated fixed-width fields directly into final BCF
payload buffers, or specialize complete row executors that combine parse,
validation, and encode rather than only replacing the opcode switch.

## Follow-Up: Direct Payload Sinks

The next pass tested direct final-buffer output for fields whose BCF
representation is known before parsing:

- exact `GT2` writes directly to a final `INT8` payload instead of scratch
  `int32_t` values plus `bcf_enc_vint()`;
- exact `AB` writes directly to a final float payload;
- strict shape executors direct-write `GT2` and optional leading `FLOAT1`
  payloads, with rollback on fallback;
- exact AD/DP/GQ/PL also carry integer range metadata into a known-range encoder
  to avoid the range pass in `bcf_enc_vint()`.

Correctness remained byte-identical for the edge fixture and 10k CCDG exact and
interpreter modes.

Parse-heavy 10k CCDG reference:

| Conversion | Baseline | Exact kernels | Direct-sink interp |
|---|---:|---:|---:|
| VCF -> uncompressed BCF | 2.51-2.68 s | 1.57-1.58 s | 2.29-2.39 s |

Full 10k CCDG compressed matrix, real seconds:

| Conversion | Baseline | Exact kernels | Direct-sink interp |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz | 9.51 | 8.53 | 9.28 |
| BCF -> BCF.gz | 7.46 | 7.46 | 7.46 |
| BCF -> VCF.gz | 11.95 | 12.00 | 12.02 |
| VCF.gz -> VCF.gz | 14.16 | 12.95 | 13.62 |
| VCF.gz -> uncompressed BCF | 2.95 | 1.92 | 2.67 |

The direct sinks are safe but small on this dataset.  The known-range encoder
was also byte-identical but did not produce a clear timing win, suggesting that
range tracking during parse still mostly trades one cost for another.  Broader
direct integer output likely needs either a cheap type-prediction/rollback
strategy or complete fused row executors that avoid both scratch traffic and
post-parse encoding for multiple fields at once.

## Follow-Up: Optimistic Guards

The fast paths now have a small circuit breaker in the cached plan state.  This
is tuned for the practical expectation that files are piecewise fixed-format,
with occasional weird rows rather than uniformly weird records.

The fast parser still validates as it parses and immediately rolls back on any
mismatch.  The new guard only decides whether to keep trying that fast parser on
later records:

- a success resets the consecutive-miss streak;
- isolated weird rows fall back once and do not disable the fast path;
- eight consecutive misses pause the fast path;
- after 128 attempts, more than 10% fallbacks also pauses it;
- paused paths cool down for 256 skipped records, then re-probe so later
  fixed-format regions can recover.

The clean CCDG path is unchanged: on the 10k subset, exact mode still reports
`10000 hits / 0 fallbacks`.  The edge fixture remains byte-identical and keeps
the expected mixed behavior:

```text
./test/test_format_plan.sh
vcf-format-plan attempts=14 hits=11 fallback=3 parsed_samples=33
vcf-format-plan attempts=14 hits=14 fallback=0 parsed_samples=42
```

The full compressed matrix was not re-recorded for this guard-only change
because the machine was under unrelated CPU load during the interrupted run.
The parse-heavy 10k BCF outputs were re-compared byte-for-byte for exact and
interpreter modes.

## Follow-Up: Generic Strict Numeric Executor

The next iteration removed the two hard-coded shape executors and replaced them
with a generic strict fixed-numeric executor.  This is the dynamic-exact version
of the FORMAT planner:

- the executor is keyed by resolved row op kinds and widths, not FORMAT field
  names;
- any fixed-width numeric op sequence is eligible;
- leading `GT2` and scalar float fields are written directly into the final BCF
  `indiv` buffer;
- integer fields carry min/max range metadata from parse into encode so
  `bcf_enc_vint()` does not rescan scratch arrays;
- any mismatch rolls back direct writes and falls back to the measured-width
  general planner or legacy parser.

Correctness checks:

```text
make -j4 test/test_view
./test/test_format_plan.sh
cmp baseline/exact/interp BCF outputs for /tmp/ccdg_chr22_10k.vcf.gz
cmp baseline/exact/interp compressed BCF outputs for /tmp/ccdg_chr22_10k.vcf.gz
```

The mixed edge fixture remains byte-identical.  It now includes reordered
numeric FORMAT fields, a scalar float away from the first FORMAT positions,
non-CCDG fixed-width numeric tags, and integer values that cross BCF int8/int16
encoding thresholds:

```text
vcf-format-plan attempts=14 hits=11 fallback=3 parsed_samples=33
vcf-format-plan attempts=14 hits=14 fallback=0 parsed_samples=42
```

Parse-heavy 10k CCDG, VCF.gz to uncompressed BCF, real seconds:

| Mode | Run 1 | Run 2 | Run 3 |
|---|---:|---:|---:|
| Baseline | 2.86 | 2.87 | 2.85 |
| Exact kernels | 1.85 | 1.85 | 1.86 |
| Dynamic strict/interp | 1.87 | 1.88 | 1.88 |

After removing the old shape-specific templates, a cleanup check still showed
exact and dynamic strict essentially tied:

| Mode | Real seconds |
|---|---:|
| Exact kernels | 1.89 |
| Dynamic strict/interp | 1.87 |

Single-run compressed VCF.gz to compressed BCF.gz, real seconds:

| Mode | Real seconds |
|---|---:|
| Baseline | 10.08 |
| Exact kernels | 9.01 |
| Dynamic strict/interp | 8.58 |

The compressed result should be read as directional because compression noise is
larger, but the outputs were byte-identical and the main parse-heavy result is
the key signal: the dynamic strict path is now within measurement noise of the
hand-written exact CCDG kernel without matching on CCDG field names.

## Follow-Up: Subagent Review Cleanup

Three review passes suggested tightening the generic executor rather than adding
more field-name cases.  The resulting cleanup:

- rolls back `v->indiv` on hard errors after direct writes, not only on shape
  fallback;
- skips range updates for padding vector-end sentinels while preserving the
  `bcf_enc_vint()` min/max contract;
- precomputes per-op base pointers and strides before the sample loop;
- removes the stale non-range encoder variant;
- expands the edge fixture with reordered and non-CCDG numeric FORMAT rows.

Post-cleanup parse-heavy 10k CCDG, VCF.gz to uncompressed BCF, real seconds:

| Mode | Run 1 | Run 2 | Run 3 |
|---|---:|---:|---:|
| Baseline | 2.83 | 2.84 | 2.86 |
| Exact kernels | 1.85 | 1.89 | 1.86 |
| Dynamic strict/interp | 1.86 | 1.83 | 1.83 |

## Follow-Up: Broader Operation Checks

After the dynamic strict executor cleanup, the 10k CCDG conversion matrix was
rerun with `test/test_view`.  Outputs were compared byte-for-byte against the
baseline output for every exact/interp cell.

Single-run format conversion matrix, real seconds:

| Conversion | Baseline | Exact kernels | Dynamic strict/interp |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz | 8.73 | 7.78 | 8.58 |
| BCF -> BCF.gz | 6.85 | 6.92 | 7.02 |
| BCF -> VCF.gz | 11.18 | 11.22 | 11.15 |
| VCF.gz -> VCF.gz | 13.26 | 12.34 | 13.01 |
| VCF.gz -> uncompressed BCF | 2.83 | 1.85 | 2.58 |

The `VCF.gz -> uncompressed BCF` interp cell above was a noisy outlier; a direct
focused rerun of that same parse-heavy case reproduced exact-speed dynamic
strict behavior:

| Mode | Run 1 | Run 2 | Run 3 |
|---|---:|---:|---:|
| Exact kernels | 1.84 | 1.82 | 1.85 |
| Dynamic strict/interp | 1.84 | 1.84 | 1.85 |

Read-only scan with `test_view -B` isolates input decode/parse without output
formatting or compression:

| Input | Mode | Run 1 | Run 2 | Run 3 |
|---|---|---:|---:|---:|
| VCF.gz | Baseline | 2.59 | 2.61 | 2.58 |
| VCF.gz | Exact kernels | 1.62 | 1.62 | 1.65 |
| VCF.gz | Dynamic strict/interp | 1.62 | 1.63 | 1.62 |
| BCF | Baseline | 0.62 | 0.61 | 0.61 |
| BCF | Exact kernels | 0.61 | 0.61 | 0.63 |
| BCF | Dynamic strict/interp | 0.61 | 0.62 | 0.61 |

Threaded compressed output with `test_view -@ 4` makes the parser win visible
again even for compressed-to-compressed workflows:

| Conversion | Baseline | Exact kernels | Dynamic strict/interp |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz, `-@ 4` | 2.64 | 2.03 | 2.06 |
| VCF.gz -> VCF.gz, `-@ 4` | 3.96 | 3.03 | 3.02 |

BCF-input conversions and BCF read-only scans remain unchanged, as expected,
because the optimization only affects VCF FORMAT parsing.

An attempted bcftools rebuild against this htslib worktree failed at link time
because the sibling bcftools checkout expects `bcf_write_take_ownership`, which
is not present in this htslib worktree.  Operation-level bcftools timings should
therefore be rerun only after pairing this htslib branch with a matching
bcftools revision or porting that API.

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
