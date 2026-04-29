# FORMAT Plan Parser Spec

This document describes the current experimental `HTS_VCF_FORMAT_PLAN` VCF
FORMAT parser and the direction for making it more general.

## Goal

Keep the existing htslib FORMAT parser as the source of truth, but add
opportunistic fast paths for repeated FORMAT layouts.  A fast path may only claim
a record when it can produce byte-identical BCF.  Otherwise it must return `-3`
and let the existing parser handle the row.

## Current Architecture

`HTS_VCF_FORMAT_PLAN=1` enables a tiered parser:

1. Handwritten exact kernels for the four dominant CCDG FORMAT layouts:
   `GT:AB:AD:DP:GQ:PL`, `GT:AD:DP:GQ:PL`,
   `GT:AB:AD:DP:GQ:PGT:PID:PL`, and `GT:AD:DP:GQ:PGT:PID:PL`.
2. A dynamic general FORMAT planner keyed by the literal FORMAT column and
   header pointer.  It resolves field IDs/types once and executes row-specific
   operations for GT, integer vectors, float vectors, and strings.
3. The existing generic htslib FORMAT parser for unsupported or suspicious
   rows.

`HTS_VCF_FORMAT_PLAN=interp` or `HTS_VCF_FORMAT_PLAN=general` skips the exact
CCDG kernels and runs only the dynamic general planner.  This mode is useful for
isolating how much performance the general approach has captured.

## Measured State

On the 10k CCDG subset, the exact tier is currently the large win.  A clean
sanity rerun on 2026-04-29 showed:

| Mode | VCF.gz read-only | VCF.gz -> uncompressed BCF |
|---|---:|---:|
| Baseline | 2.58 s | 2.83 s |
| `HTS_VCF_FORMAT_PLAN=1` | 1.61 s | 1.86 s |
| `HTS_VCF_FORMAT_PLAN=interp` | 2.34 s | 2.55 s |

The earlier docs overstated the dynamic strict/interpreter result.  The dynamic
planner is correct and modestly faster than baseline, but it does not yet match
the handwritten CCDG kernels.

The next development target is to move the exact-kernel advantages into a
dynamic shape executor so common fixed-format regions can get exact-like speed
without field-name-specific kernels.

## Correctness Rules

The planned parser must preserve these invariants:

- No planned parsing while `h->keep_samples` is active.
- Header IDs and types are resolved before execution.
- Duplicate FORMAT tags use the generic parser.
- Undefined tags use the generic parser, preserving current dummy-header
  behavior and warnings.
- GT encoding must match generic htslib phasing semantics, including haploid
  genotypes, missing alleles, multidigit allele indexes, and VCF 4.4 prefix
  phasing.
- Numeric vectors use observed or provably fixed row width and pad shorter
  samples with vector-end sentinels.
- Strings use observed maximum byte length and zero-pad shorter samples.
- Integer and float overflow/error behavior must either match generic htslib or
  force fallback.
- Any fast path that writes directly into `v->indiv` must save the original
  length and roll back before fallback.

## Dynamic Planner

The general planner compiles the literal FORMAT string into a cached op list.
After seeing a record, it resolves the ops to row-local opcodes such as `GT2`,
`GT`, `INT1`, `INT2`, `INT3`, `INTN`, `FLOAT1`, `FLOATN`, and `STR`.

For rows whose widths can be predicted from the header and allele count, the
planner first tries a strict numeric executor.  That path validates shape while
parsing, carries integer min/max metadata into BCF integer encoding, and can
direct-write a leading `GT2`/`FLOAT1` prefix.  If the row is sparse, stringy,
malformed, or otherwise not byte-identical, it falls back to the measured-width
general planner.

Today, the strict/general path still has enough overhead that it trails the
handwritten CCDG kernels on the CCDG benchmark.  Likely remaining gaps include
per-field dispatch, measured-width/string handling for `PGT/PID`, scratch-buffer
traffic, and generic encode costs.

## Guard Policy

Each cached exact/general plan has a small runtime guard:

- attempts, hits, fallbacks,
- consecutive miss streak,
- temporary cooldown.

An isolated fallback does not disable the fast path.  A plan is paused only
after eight consecutive misses, or after at least 128 attempts with more than
10% fallbacks.  After 256 skipped records, the plan probes again so later
fixed-format regions can recover the optimized path.

For exact CCDG kernels, a paused exact guard routes the row to the dynamic
general planner.  For general plans, a paused strict guard skips directly to the
measured-width planner, and a paused general guard returns to legacy htslib
parsing.

## Edge Fixture

`test/format-plan-edge.vcf` is CCDG-shaped but includes awkward realistic rows:

- the exact CCDG layouts,
- reordered FORMAT fields,
- non-CCDG numeric tag names with fixed widths,
- integer values around BCF int8/int16 boundaries,
- multiallelic AD/PL and GL,
- haploid GT,
- multidigit allele indexes,
- fixed integer vectors,
- string FORMAT fields,
- exact-kernel fallbacks that the dynamic planner can still handle.

Run:

```sh
./test/test_format_plan.sh
```

The script writes BCF through the generic parser, `HTS_VCF_FORMAT_PLAN=1`, and
`HTS_VCF_FORMAT_PLAN=interp`, then compares the outputs with `cmp`.

## Next Work

- Make a dynamic fixed-shape executor that captures the CCDG exact-kernel wins
  without matching on field names.
- Specialize common string-bearing shapes such as `PGT/PID` without baking in
  CCDG tag names.
- Reduce per-sample opcode dispatch in hot FORMAT shapes.
- Expand direct final-buffer output only where BCF type selection remains
  byte-identical or can cheaply roll back.
- Keep the exact kernels as a performance oracle while iterating, then remove
  or demote them once the dynamic executor catches up.
