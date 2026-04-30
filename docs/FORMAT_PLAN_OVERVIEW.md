# Dynamic FORMAT Plan Overview

This branch adds an optional fast path for parsing VCF `FORMAT` sample columns.
The goal is to speed up common, header-described FORMAT layouts without writing
one-off kernels for exact FORMAT strings such as `GT:AD:DP:GQ:PL`.

## What It Does

When `HTS_VCF_FORMAT_PLAN` is enabled, htslib first tries to compile the record's
literal FORMAT string into a small list of per-tag operations.  The plan is
driven by the active VCF header: each tag contributes its key, type, declared
number model, and whether the current row needs width measurement.

Compiled plans live in private header-owned cache state.  The cache is cleared
when the header dictionaries are resynchronised, and the optimized path declines
to run while the header has unsynced mutations.  That keeps cached supported and
unsupported decisions tied to the exact header metadata that produced them.

If the row fits the supported operation set, the dynamic executor parses samples
and writes BCF's transposed FORMAT layout directly.  If anything looks unsafe or
unsupported, htslib falls back to the generic parser for the whole FORMAT
column.  The planner also keeps a small profitability gate: schemas dominated by
measured strings plus float vectors, such as `GT:FT:PID:GL:DP`, currently use
the generic parser because the dynamic path's width-measurement work costs
more than it saves.

The optimized path also supports selected-sample reads.  When
`bcf_hdr_set_samples()` is active, it scans the original sample columns, skips
unretained samples, and writes the retained samples densely into the BCF FORMAT
blocks.

Fallbacks are whole-row, but they are now classified for diagnostics when
`HTS_VCF_FORMAT_PLAN_STATS=1` is set.  The current reason counters distinguish
unsupported schemas, guard cooldowns, numeric-width limits, string-width limits,
GT shape misses, parse failures, separator mismatches, and sample-count
mismatches.

## Why This Shape

The important design choice is tag-level composition.  A file does not need an
exact hardcoded FORMAT string to benefit.  For example, these can all share the
same dynamic machinery when their tags are described by supported header
metadata:

- `GT:AD`
- `GT:AD:DP:PL`
- `GT:AB:AD:DP:GQ:PGT:PID:PL`
- reordered numeric/string tags
- supersets with additional supported tags

This is deliberately more general than the earlier experimental exact kernels.
Those kernels were fast, but brittle: adding or removing one tag could miss the
optimized path entirely.

## Where It Helps

The feature is most useful for sample-rich VCF text input where FORMAT parsing is
a meaningful part of total runtime:

- large `GT`-only genotype VCFs;
- likelihood-heavy VCFs with fields such as `AD`, `PL`, `DP`, `GQ`, `AB`, and
  phase strings;
- conversion paths such as VCF.gz to BCF where text FORMAT parsing is exposed;
- workloads with repeated FORMAT layouts across many records.

In the latest bcftools-style timing, the real 1000G chr22 GT workload sped up
from 27.48 s to 8.99 s unthreaded, and from 26.71 s to 6.94 s at 4 threads.
The likelihood-heavy synthetic workload improved more modestly, from 4.43 s to
3.94 s unthreaded and from 3.47 s to 3.02 s at 4 threads.

With bcftools selecting the first two samples from each input, the same 1000G
GT workload improved from 26.51 s to 9.77 s unthreaded and from 25.99 s to
8.84 s at 4 threads.  Selected-sample likelihood-heavy rows are still faster,
but the gains are smaller because much less FORMAT payload is emitted.

Broader bcftools commands follow the same pattern.  `bcftools view`,
`bcftools query` of FORMAT values, and genotype filters benefit when they expose
sample FORMAT parsing.  Site-only queries, `stats`, and `merge` are mostly
neutral because their runtime is dominated by non-FORMAT work, output writing,
or command-level bookkeeping.  A controlled `bcftools merge` self-merge check
produced byte-identical output and was neutral-to-positive across the small
merge manifest.

## Drawbacks

The MVP intentionally keeps fallback whole-row.  It does not parse supported
tags dynamically while delegating only one unsupported tag to the generic
parser.  That makes correctness easier to reason about, but a single unsupported
tag or malformed row means the entire FORMAT column uses the generic parser.

Known fallback cases include:

- undefined FORMAT tags that require production header repair;
- unsupported header types or number models;
- low-profit string/float-heavy schemas;
- duplicate FORMAT tags;
- malformed separators or unexpected sample cardinality;
- row-local widths above the bounded fast-path limit;
- GT encodings outside the simple fast-path representation.

The path is also not always faster.  Some string/float-heavy layouts are roughly
at parity or slightly slower than baseline because the dynamic path still pays
measurement, dispatch, and scratch-buffer costs.

The current planned width limits are intentionally conservative: measured
numeric vectors are capped at 64 values, and measured strings are capped at
256 bytes.  Rows above those limits use the generic parser; numeric/string width
misses do not by themselves disable the schema for later rows.

Correctness checks for this path now live in the normal htslib test harness, not
only in the benchmark directory.  `make check` runs black-box byte-identity
fixtures through `test/test.pl`, selected-sample checks, malformed-input checks,
and focused header-cache generation coverage.

## User-Facing Controls

```text
unset / 0       generic parser only
1              dynamic per-tag planner, then generic fallback
```

The benchmark harness reports only `HTS_VCF_FORMAT_PLAN=1` as `plan`.
Other values are treated as disabled.

## Related Docs

- `docs/FORMAT_PLAN_CURRENT.md`: current implementation, supported shapes,
  correctness rules, and benchmark tables.
- `docs/FORMAT_PLAN_EXPERIMENT_LOG.md`: chronological log of approaches tried,
  results, reversions, and retained lessons.
