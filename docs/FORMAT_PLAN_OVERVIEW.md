# VCF FORMAT Planner Review Notes

This is an interim review note for the opt-in VCF FORMAT planner in `vcf.c`.
It is intended to make the implementation easier to review in this branch; it
is not proposed as permanent user-facing documentation.

## Purpose

The planner is an optional fast path for parsing VCF sample `FORMAT` columns
into BCF.  It is disabled by default and only runs when:

```sh
HTS_VCF_FORMAT_PLAN=1
```

Unset, `0`, or unknown values use the existing generic parser.  The hard
correctness rule is byte-identical BCF output compared with the generic parser.
Any unsupported or suspicious FORMAT column falls back for the whole column.

`HTS_VCF_FORMAT_PLAN_STATS=1` enables temporary diagnostic counters for review
and benchmarking.  These environment variables are branch controls, not stable
public API.

## Entry Point

`vcf_parse_format()` tries `vcf_parse_format_planned()` before the generic
FORMAT parser when the environment gate is enabled.  The planned parser returns
success only after it has fully emitted the FORMAT column.  It returns the
fallback code before the caller invokes the generic parser when compilation,
width resolution, parsing, or row validation does not meet the supported
contract.

The planned path never commits partial FORMAT output after detecting a row-local
failure.  Rollback and fallback are part of the normal control flow.

## Plan Cache

Plans are stored in private `bcf_hdr_aux_t` state.  Cache keys are the literal
FORMAT string plus the active private header generation.  The cache stores both
supported and unsupported plans so repeated unsupported schemas avoid repeated
tokenization and metadata lookup.

`bcf_hdr_sync()` clears cached plans and advances the generation because FORMAT
ids, types, and lengths are header-local.  The planner refuses dirty headers.

## Compilation

Compilation works from header metadata rather than exact whole-FORMAT string
kernels.  For each FORMAT token, the compiler records:

- the header id;
- declared type;
- declared number model;
- whether the row needs record-local width resolution or sample-text scanning;
- the executor row kind.

The compiler rejects empty tokens, undefined tags, duplicate tags, unsupported
types or number models, and non-standard `GT` declarations.  Tokenization does
not collapse empty fields, so malformed schemas such as `GT::DP` still fall
back in a way that preserves generic-parser behavior.

## Supported Rows

The current executor has six row kinds:

- `GT2`
- `INT1`
- `INTVEC`
- `FLOAT1`
- `FLOATN`
- `STR`

The earlier width-specific integer-vector row kinds were removed.  A single
`INTVEC` path handles fixed and row-local integer widths, including the
over-width comma check needed to preserve fallback behavior.

Supported shapes include:

- simple diploid `GT` values with one-character alleles or missing values,
  separated by `/` or `|`, including phased-missing forms such as `.|.`,
  `0|.`, and `.|0`;
- integer and float scalar fields;
- integer and float vector fields within the planner width cap;
- numeric `Number=A`, `Number=R`, and `Number=G` widths resolved from the
  current record allele count;
- bounded measured `Number=.` numeric rows;
- bounded `Type=String,Number=1` rows;
- selected-sample parsing via `bcf_hdr_set_samples()`.

Unsupported or intentionally generic cases include undefined tags, duplicate
FORMAT tags, dirty headers, unsupported type or number declarations, unsupported
GT encodings, malformed separators, unsafe row widths, and string/float-heavy
layouts that do not benefit from planning.

## Width Resolution

Header-fixed rows use the declared width directly, after resolving
allele-dependent widths for the current record.  Numeric widths must fit the
planner cap of 64 values.

Measured numeric and string rows perform a first pass over original sample
columns.  This is required to match the generic parser's width and padding
rules.  If samples are selected, the planner still scans original sample
columns but measures and emits retained samples densely.

Strings are capped at 256 bytes in the planned path.  Wider string rows fall
back for the whole FORMAT column.

## Fallback Contract

Fallback is expected and intentional.  It happens before generic parsing when
the compiler or row executor sees unsupported structure, unsupported widths,
unexpected separators, unsupported GT shape, parse failures, sample-count
mismatch, or allocation/internal consistency errors.

Diagnostic fallback reasons are:

- `unsupported`
- `numeric_width`
- `string_width`
- `gt_shape`
- `parse`
- `separator`
- `sample_count`

## Test And Benchmark Evidence

Focused correctness checks used for this review branch:

```sh
make test/test_view test/test_format_plan_cache bgzip tabix
./test/test_format_plan_cache
perl test/test.pl -F test_vcf_format_plan
test/maintainer/check_spaces.pl vcf.c docs/FORMAT_PLAN_OVERVIEW.md \
  test/format-plan-malformed-fields.vcf test/test.pl
git diff --check
```

The focused FORMAT-plan test fragment covers disabled/unknown environment
behavior, selected samples, malformed FORMAT tokens, malformed numeric fields,
phased missing GT values, cache invalidation after header metadata changes, and
fallback after partial planned parsing.

The current public-fork PR body contains the maintainer-facing performance
summary.  The compact benchmark artifacts live on the corpus branch
`feature/vcf-parsing-speedup-corpus`, including the current `test_view` and
bcftools summaries for commit `bd643182c8fa722abbc0cb89860263a90bb97020`.
