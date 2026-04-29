# FORMAT Plan Parser Spec

This document describes the current `HTS_VCF_FORMAT_PLAN` VCF FORMAT parser.
The older exact CCDG kernels and dynamic likelihood-shape tier were removed;
the optimized path is now a single dynamic per-tag planner with production
fallback.

## Goal

Keep the existing htslib FORMAT parser as the source of truth, while adding an
opportunistic fast path for repeated, header-described FORMAT layouts.  The fast
path may only claim a row when it can produce byte-identical BCF.  Otherwise it
must return `-3` and let the existing parser handle the whole FORMAT column.

## Current Architecture

`HTS_VCF_FORMAT_PLAN` controls the planned parser:

```text
unset / 0       production parser only
1              dynamic per-tag planner, then production fallback
interp/general same dynamic per-tag planner, then production fallback
```

All enabled spellings now use the same implementation.  The benchmark harness
may still run both `1` and `interp` as separate modes, but they are intended to
match except for normal timing noise.

The planned parser has four stages:

1. Compile the literal FORMAT string and active header into a cached list of
   per-tag operations.
2. Resolve row-local widths from header `Number` metadata, allele count, and a
   bounded measurement pass for strings or `Number=.` numeric vectors.
3. Parse sample fields into BCF's transposed FORMAT layout with a composable
   executor.
4. Fall back to the production FORMAT parser for unsupported or suspicious rows.

## Supported Tags

The planner is tag-composable rather than full-string-specialized.  It can claim
layouts such as `GT:AD`, `GT:AD:DP:PL`,
`GT:AB:AD:DP:GQ:PGT:PID:PL`, reordered fields, and supersets when each tag has
supported header metadata.

Supported FORMAT tag shapes:

- `GT` declared as `Type=String,Number=1`, with simple diploid encodings on the
  fast path.
- Integer fields with fixed `Number=N`, `Number=A`, `Number=R`, `Number=G`, or
  bounded measured `Number=.` row widths.
- Float fields with the same number models as integer fields.
- String fields declared as `Type=String,Number=1`, measured per row.

Unsupported tags or unsupported row-local encodings fall back whole-row.

## Correctness Rules

The planned parser must preserve these invariants:

- No planned parsing while `h->keep_samples` is active.
- Header IDs, types, and number models are resolved before execution.
- Duplicate FORMAT tags use the production parser.
- Undefined tags use the production parser, preserving dummy-header behavior and
  warnings.
- GT encoding must match htslib phasing semantics; encodings outside the simple
  fast path must force fallback.
- Numeric vectors use observed or provably fixed row width and pad shorter
  samples with vector-end sentinels.
- Strings use observed maximum byte length and zero-pad shorter samples.
- Integer and float overflow/error behavior must either match production htslib
  or force fallback.
- Any fast path that writes directly into `v->indiv` must save the original
  length and roll back before fallback.

## Dynamic Planner

The planner compiles the literal FORMAT string into cached opcodes keyed by
header pointer plus FORMAT text.  Header-local ids and type metadata make plans
unsafe to share across headers.

After seeing a record, it resolves the reusable op list to row-local operations
such as `GT2`, `INT1`, `INT2`, `INT3`, `INTN`, `FLOAT1`, `FLOATN`, and `STR`.
`Number=A`, `Number=R`, and `Number=G` widths come from the current allele
count.  String and `Number=.` numeric widths are measured across the row before
execution.

The executor writes BCF's transposed FORMAT layout.  Leading fixed-width
`GT2`/`FLOAT1` rows can be written directly into `v->indiv`; other rows are
staged in header scratch memory and encoded after sample parsing so integer
range and observed-width metadata are known.

## Guard Policy

Each cached dynamic plan has a small runtime guard:

- attempts, hits, fallbacks,
- consecutive miss streak,
- temporary cooldown.

An isolated fallback does not disable the fast path.  A plan is paused only
after eight consecutive misses, or after at least 128 attempts with more than
10% fallbacks.  After 256 skipped records, the plan probes again so later
fixed-format regions can recover the optimized path.

## Tests

`./test/test_format_plan.sh` writes BCF through:

- the production parser,
- `HTS_VCF_FORMAT_PLAN=1`,
- `HTS_VCF_FORMAT_PLAN=interp`.

It compares the planned outputs against baseline with `cmp`.  The fixtures cover
subsets, supersets, reordered fields, measured numeric fields, strings,
malformed header shapes, and deliberate row-local fallback cases.

## Next Work

- Add selected-sample support so `keep_samples` does not require whole-row
  fallback.
- Reduce per-sample opcode dispatch in hot FORMAT layouts.
- Improve string and measured-width handling without losing byte identity.
- Consider a later executor-generation layer if generic per-op dispatch remains
  the main gap to historical exact-kernel speed.
