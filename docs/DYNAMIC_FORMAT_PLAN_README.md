# Dynamic FORMAT Plan

This branch adds an optional dynamic fast path for parsing VCF `FORMAT` sample
columns.  The goal is to speed up common, header-described FORMAT layouts
without hardcoding exact full FORMAT strings such as `GT:AD:DP:GQ:PL`.

## How It Works

When `HTS_VCF_FORMAT_PLAN` is enabled, `vcf_parse_format()` first tries to
compile the record's literal FORMAT string into a small list of per-tag
operations.  Compilation uses the active header, so the plan records each tag's
header key, type, declared number, and whether the row needs width measurement.

The executor then parses samples with that op list and writes BCF's transposed
FORMAT layout.  If compilation or row-local validation fails, it returns to the
existing production parser for the whole FORMAT column.

Supported environment values:

- unset or `0`: use the production parser only.
- `1`, `interp`, or `general`: use the dynamic plan, with production fallback.

The old exact FORMAT kernels and optional SIMD tab-scanning front-end have been
removed.  All enabled spellings now route through the same dynamic path.

## Supported Cases

The fast path is tag-composable rather than full-string-specialized.  It can
handle subsets, reordered fields, and supersets when each tag is described by
header metadata that the executor supports.

Currently supported FORMAT tag shapes:

- `GT` declared as `Type=String,Number=1`, with simple diploid encodings on the
  fast path.
- Integer fields with fixed `Number=N`, `Number=A`, `Number=R`, `Number=G`, or
  bounded measured `Number=.` row widths.
- Float fields with the same number models as integer fields.
- String fields declared as `Type=String,Number=1`, measured per row.

Examples that can use the dynamic path include `GT:AD`, `GT:AD:DP:PL`,
`GT:AB:AD:DP:GQ:PGT:PID:PL`, and reordered/superset layouts with additional
supported tags.

## Fallback Behavior

Fallback is intentionally whole-row for the MVP.  The dynamic parser does not
mix optimized handling for some tags with production handling for other tags in
the same FORMAT column.  This keeps BCF layout, warning behavior, and error
recovery aligned with the existing parser when a row is unusual.

Known fallback cases include:

- sample subsetting via `keep_samples`;
- undefined FORMAT tags that require production header repair;
- unsupported header types or number models;
- malformed sample separators or unexpected sample cardinality;
- row-local widths above the current bounded fast-path limit;
- GT encodings outside the simple fast-path representation.

## Tests And Benchmarks

Focused validation lives in `test/test_format_plan.sh`, including byte-for-byte
comparisons between production parsing and planned parsing on edge-case FORMAT
fixtures.

The larger benchmark corpus lives under `bench/format-shape/large/`.  The
benchmark script runs `baseline` and `plan` modes.  `plan` is
`HTS_VCF_FORMAT_PLAN=1`; `interp` and `general` remain accepted aliases, but
they use the same dynamic executor and are not emitted as separate timing rows.

Threaded scaling checks use `bench/format-shape/scripts/run_thread_bench.sh`.
That runner exercises representative large inputs with unthreaded execution plus
`test/test_view -@ 2`, `-@ 4`, and `-@ 8`, and writes a `threads` column in its
timing output.

Production-style checks use `bench/format-shape/scripts/run_bcftools_bench.sh`
with a bcftools binary built against this htslib tree.  It runs
`bcftools view --no-version -Ob -l 0` in baseline and planned modes over the same
representative threaded manifest.

Latest documented results are in
`docs/DYNAMIC_FORMAT_SHAPE_EXECUTOR_SCRATCHPAD.md`.
