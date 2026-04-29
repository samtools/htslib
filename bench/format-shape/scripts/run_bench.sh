#!/bin/sh
set -eu

test_view=${TEST_VIEW:-./test/test_view}
inputs=${1:-bench/format-shape/inputs.tsv}
outdir=${OUTDIR:-bench/format-shape/results}
mkdir -p "$outdir"

timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"

printf 'name\tmode\treal\tuser\tsys\tattempts\thits\tfallback\tparsed_samples\tshape_attempts\tshape_hits\tshape_fallback\n' > "$timings"
printf 'name\tcomparison\tstatus\n' > "$checks"

tail -n +2 "$inputs" | while IFS='	' read -r name path source
do
    base_out="$outdir/$name.baseline.bcf"
    exact_out="$outdir/$name.exact.bcf"
    interp_out="$outdir/$name.interp.bcf"

    for mode in baseline exact interp
    do
        err="$outdir/$name.$mode.stderr"
        out="$outdir/$name.$mode.bcf"
        case "$mode" in
            baseline)
                env HTS_VCF_FORMAT_PLAN=0 /usr/bin/time -p "$test_view" -b -l 0 "$path" > "$out" 2> "$err"
                ;;
            exact)
                env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 HTS_VCF_FORMAT_PLAN_SHAPE_STATS=1 \
                    /usr/bin/time -p "$test_view" -b -l 0 "$path" > "$out" 2> "$err"
                ;;
            interp)
                env HTS_VCF_FORMAT_PLAN=interp HTS_VCF_FORMAT_PLAN_STATS=1 HTS_VCF_FORMAT_PLAN_SHAPE_STATS=1 \
                    /usr/bin/time -p "$test_view" -b -l 0 "$path" > "$out" 2> "$err"
                ;;
        esac

        awk -v name="$name" -v mode="$mode" '
            /^real / { real=$2 }
            /^user / { user=$2 }
            /^sys / { sys=$2 }
            /^vcf-format-plan / {
                for (i=1; i<=NF; i++) {
                    split($i, kv, "=")
                    if (kv[1] == "attempts") attempts=kv[2]
                    else if (kv[1] == "hits") hits=kv[2]
                    else if (kv[1] == "fallback") fallback=kv[2]
                    else if (kv[1] == "parsed_samples") parsed=kv[2]
                }
            }
            /^vcf-format-likelihood-shape / {
                for (i=1; i<=NF; i++) {
                    split($i, kv, "=")
                    if (kv[1] == "attempts") shape_attempts=kv[2]
                    else if (kv[1] == "hits") shape_hits=kv[2]
                    else if (kv[1] == "fallback") shape_fallback=kv[2]
                }
            }
            END {
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                       name, mode, real+0, user+0, sys+0,
                       attempts+0, hits+0, fallback+0, parsed+0,
                       shape_attempts+0, shape_hits+0, shape_fallback+0
            }
        ' "$err" >> "$timings"
    done

    if cmp "$base_out" "$exact_out" >/dev/null 2>&1; then
        printf '%s\tbaseline_vs_exact\tok\n' "$name" >> "$checks"
    else
        printf '%s\tbaseline_vs_exact\tDIFF\n' "$name" >> "$checks"
    fi
    if cmp "$base_out" "$interp_out" >/dev/null 2>&1; then
        printf '%s\tbaseline_vs_interp\tok\n' "$name" >> "$checks"
    else
        printf '%s\tbaseline_vs_interp\tDIFF\n' "$name" >> "$checks"
    fi
done

printf 'wrote %s and %s\n' "$timings" "$checks"
