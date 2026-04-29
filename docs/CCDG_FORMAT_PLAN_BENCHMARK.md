# CCDG FORMAT Plan Benchmark Checkpoint

Date: 2026-04-29

Worktree: `/tmp/htslib-vcf-avx-sanity`

Branch: `codex/vcf-avx-sanity`

## Current Takeaway

The experimental FORMAT planner is viable, but the current large CCDG win comes
from the handwritten exact CCDG kernels, not yet from the fully dynamic
strict/interpreter path.

The dynamic general planner is correct and modestly faster than baseline.  It is
the path we want to improve next, using the exact kernels as a performance
oracle.

## Modes

```sh
HTS_VCF_FORMAT_PLAN=0       # baseline generic parser
HTS_VCF_FORMAT_PLAN=1       # exact CCDG kernels, then dynamic general fallback
HTS_VCF_FORMAT_PLAN=interp  # dynamic general planner only
HTS_VCF_FORMAT_PLAN_STATS=1 # print planner counters from test/test_view
```

## Data

Source file:

```text
/Users/jeremiah.li/geneticoptims/inplace-htslib-refactor/data/original/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

Subset used for the current benchmark:

```text
/tmp/ccdg_chr22_10k.vcf.gz
/tmp/ccdg_chr22_10k.bcf
```

The subset contains 10,000 variant records and 3,202 samples.  The observed
FORMAT distribution is:

| Records | FORMAT |
|---:|---|
| 4,681 | `GT:AB:AD:DP:GQ:PL` |
| 3,774 | `GT:AB:AD:DP:GQ:PGT:PID:PL` |
| 813 | `GT:AD:DP:GQ:PL` |
| 732 | `GT:AD:DP:GQ:PGT:PID:PL` |

The exact CCDG tier covers all four layouts.

## Clean Sanity Rerun

These numbers were rerun after noticing that an earlier table mislabeled the
dynamic/interpreter result.  Timings are single wall-clock runs on the 10k CCDG
subset, so treat them as directional.

| Mode | VCF.gz read-only | VCF.gz -> uncompressed BCF |
|---|---:|---:|
| Baseline | 2.58 s | 2.83 s |
| Exact + dynamic fallback | 1.61 s | 1.86 s |
| Dynamic general only | 2.34 s | 2.55 s |

Planner counters on VCF.gz -> uncompressed BCF:

| Mode | Attempts | Hits | Fallback | Parsed samples |
|---|---:|---:|---:|---:|
| Exact + dynamic fallback | 10,000 | 10,000 | 0 | 32,020,000 |
| Dynamic general only | 10,000 | 10,000 | 0 | 32,020,000 |

Both planned modes are byte-identical against baseline in the sanity tests, but
the exact tier is much faster.

## Broader Conversion Matrix

Earlier single-run compressed conversion checks used `test/test_view` and
compared outputs byte-for-byte with `cmp`.

| Conversion | Baseline | Exact + dynamic fallback | Dynamic general only |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz | 8.73 s | 7.78 s | 8.58 s |
| BCF -> BCF.gz | 6.85 s | 6.92 s | 7.02 s |
| BCF -> VCF.gz | 11.18 s | 11.22 s | 11.15 s |
| VCF.gz -> VCF.gz | 13.26 s | 12.34 s | 13.01 s |
| VCF.gz -> uncompressed BCF | 2.83 s | 1.85 s | 2.58 s |

BCF-input conversions are unchanged, as expected, because this optimization only
affects VCF text FORMAT parsing.

Threaded compressed output with `test_view -@ 4` makes the parser win visible
even for compressed-to-compressed workflows:

| Conversion | Baseline | Exact + dynamic fallback | Dynamic general only |
|---|---:|---:|---:|
| VCF.gz -> BCF.gz, `-@ 4` | 2.64 s | 2.03 s | 2.06 s |
| VCF.gz -> VCF.gz, `-@ 4` | 3.96 s | 3.03 s | 3.02 s |

The threaded dynamic-only numbers should be rerun before drawing strong
conclusions; the clean single-thread rerun shows dynamic-only is not yet at
exact-kernel speed.

## Edge Fixture

`./test/test_format_plan.sh` compares baseline, `HTS_VCF_FORMAT_PLAN=1`, and
`HTS_VCF_FORMAT_PLAN=interp` on `test/format-plan-edge.vcf`.

Current output:

```text
vcf-format-plan attempts=14 hits=11 fallback=3 parsed_samples=33
vcf-format-plan attempts=14 hits=14 fallback=0 parsed_samples=42
```

The first line is `HTS_VCF_FORMAT_PLAN=1`: exact kernels claim the CCDG-shaped
rows and intentionally fall back for rows outside their narrow shape.  The
second line is dynamic-only: the general planner handles all 14 fixture rows.

## Profiling Notes

After `PGT:PID` support, the generic FORMAT fallback is no longer a meaningful
cost for the CCDG benchmark when exact kernels are enabled.  A macOS `sample`
profile of VCF.gz -> uncompressed BCF on the 100k subset showed the next hot
areas inside the planned path:

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

This is statistical sampling, not exact cycle accounting.  Directionally, the
next parser-side targets are integer-vector parsing, `PGT/PID` string handling,
per-sample dispatch, and repeated BCF integer encoding work.

## Checkpoint Recommendation

Commit this state as an honest experimental checkpoint:

- keep the exact CCDG kernels because they establish the upper-bound target;
- keep the dynamic general planner and edge fixture because they are the path to
  a general solution;
- keep benchmark docs explicit that dynamic-only is not yet the big win;
- do not open an upstream-facing PR until the dynamic executor closes more of
  the gap or the PR is framed as an experimental CCDG-specialized prototype.

## Next Work

The highest-value next step is to make a dynamic fixed-shape executor that
captures the exact-kernel benefits without matching on CCDG field names.  The
target is exact-like speed for piecewise fixed FORMAT regions with quick
fallback when a row leaves the proven shape.

An attempted bcftools rebuild against this htslib worktree failed at link time
because the sibling bcftools checkout expects `bcf_write_take_ownership`, which
is not present in this htslib worktree.  Operation-level bcftools timings should
be rerun only after pairing this htslib branch with a matching bcftools revision
or porting that API.
