# FORMAT Planner Current State

This document describes the current `HTS_VCF_FORMAT_PLAN=1` implementation
after simplifying the integer-vector executor.

## Entry Point

`vcf_parse_format()` calls `vcf_parse_format_planned()` before the generic
FORMAT parser.  The planned path returns:

- `0` after a successful planned parse.
- `-3` when the caller should run the generic parser.
- a hard error for allocation or internal consistency failures.

The environment gate is exact: only `HTS_VCF_FORMAT_PLAN=1` enables the planned
path.  `HTS_VCF_FORMAT_PLAN_STATS=1` records attempts, hits, fallback counts,
and parsed sample counts.

## Plan Cache

Plans are cached per header in private `bcf_hdr_aux_t` state.  Cache keys are
the literal FORMAT string plus the current private header generation.  Both
supported and unsupported FORMAT strings are cached, so repeated fallback cases
avoid recompilation.  The cache grows from 16 to 128 entries and then evicts
entries in a simple rotating order, preferring unsupported plans when possible.

`bcf_hdr_sync()` clears the cache and increments the generation because FORMAT
key ids, types, and lengths are header-local.

## Compilation

The compiler tokenizes the FORMAT string without collapsing empty tokens, so
malformed strings like `GT::DP` still fall back exactly as the generic parser
expects.  It rejects:

- empty FORMAT tokens;
- unknown or undefined tags;
- duplicate tags;
- unsupported header types;
- unsupported length models;
- non-standard `GT` declarations.

Supported non-GT types are integer, float, and string tags with enough header
metadata to reproduce the generic BCF layout.  The support is still composition
aware: measured-string plus float-vector layouts are kept on the generic parser
unless the row also contains integer-vector work.

## Row Kinds

The current row executor uses six row kinds:

- `VCF_FORMAT_ROW_GT2`
- `VCF_FORMAT_ROW_INT1`
- `VCF_FORMAT_ROW_INTVEC`
- `VCF_FORMAT_ROW_FLOAT1`
- `VCF_FORMAT_ROW_FLOATN`
- `VCF_FORMAT_ROW_STR`

The earlier width-specific integer row kinds were intentionally removed.  A
single integer-vector parser now handles fixed and row-local integer widths,
including the over-width comma check needed to preserve fallback behavior.

## Width Resolution

Header-fixed fields use the declared width as the initial row width, including
`Number=A`, `Number=R`, and `Number=G` after applying the current record allele
count.  Numeric widths must fit the planner cap of 64 values.

Measured fields perform a first pass over original sample columns.  This is
required for `Number=.` numeric rows and string spans, and it is also where the
planner validates per-sample separators.  If samples have been selected with
`bcf_hdr_set_samples()`, the planner scans original columns but only measures
and emits retained samples.

Numeric measured widths are capped at 64.  String widths are capped at 256.
Integer and float vector rows may then be compacted to the observed row maximum
when that is what the generic parser would encode.

## Fallback Contract

Any row-level parse failure falls back for the whole FORMAT column.  The
planned path must not partially keep planned output after detecting malformed
text, unsupported widths, unexpected separators, unsupported GT shape, or
sample-count mismatch.

Diagnostic fallback reasons are:

- `unsupported`
- `numeric_width`
- `string_width`
- `gt_shape`
- `parse`
- `separator`
- `sample_count`

## Tests

Focused validation commands used for the current state:

```sh
make test/test_view test/test_format_plan_cache bgzip tabix
./test/test_format_plan_cache
perl test/test.pl -F test_vcf_format_plan
git diff --check
```

The `test_vcf_format_plan` fragment compares generic, planned, disabled, and
unknown environment behavior, including selected samples and expected fallback
statistics.  It passed 21/21 tests after the simplification.

The standalone `test/test_format_plan.sh` referenced in older branch notes is
not present in this checkout; the Perl test fragment is the active focused
correctness surface here.

## Benchmark Snapshot

The local untracked `bench/format-shape/**/*.vcf.gz` corpus snapshot was run
against the pre-simplification planner at
`0999351c040793ae9763b2f9ba7550b45fd9bd11` and the simplified branch.  Both
used:

```sh
HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 test/test_view -b -l 0 INPUT >/dev/null
```

Each input had one warmup per variant followed by three paired measured runs.
All 27 streaming BCF comparisons were byte-identical.

Aggregate wall-clock delta:

| Set | n | Faster | Slower | Neutral | Delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| All inputs | 27 | 18 | 3 | 6 | -3.11% |
| Inputs with mean real >= 0.1s | 19 | 17 | 1 | 1 | -3.13% |
| FORMAT-hit inputs with mean real >= 0.1s | 15 | 14 | 0 | 1 | -3.48% |

The slower cases were tiny no-FORMAT or site-only workloads with millisecond
absolute deltas.  Among nontrivial FORMAT-hit workloads, the simplified
executor had no measured slowdown in this run.
