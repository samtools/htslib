#!/bin/sh
set -eu

bcftools=${BCFTOOLS:-bcftools}
inputs=${1:-bench/format-shape/large/threaded-inputs.tsv}
outdir=${OUTDIR:-bench/format-shape/large/results-bcftools}
keep_outputs=${KEEP_OUTPUTS:-1}
threads_list=${THREADS_LIST:-0 2 4 8}
mkdir -p "$outdir"

timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"

printf 'name\tthreads\tmode\treal\tuser\tsys\n' > "$timings"
printf 'name\tthreads\tcomparison\tstatus\n' > "$checks"

tail -n +2 "$inputs" | while IFS='	' read -r name path source
do
    for threads in $threads_list
    do
        base_out="$outdir/$name.t$threads.baseline.bcf"
        plan_out="$outdir/$name.t$threads.plan.bcf"
        thread_args=
        if [ "$threads" != 0 ]; then
            thread_args="--threads $threads"
        fi

        for mode in baseline plan
        do
            err="$outdir/$name.t$threads.$mode.stderr"
            out="$outdir/$name.t$threads.$mode.bcf"
            case "$mode" in
                baseline)
                    env HTS_VCF_FORMAT_PLAN=0 /usr/bin/time -p "$bcftools" view --no-version -Ob -l 0 $thread_args -o "$out" "$path" 2> "$err"
                    ;;
                plan)
                    env HTS_VCF_FORMAT_PLAN=1 /usr/bin/time -p "$bcftools" view --no-version -Ob -l 0 $thread_args -o "$out" "$path" 2> "$err"
                    ;;
            esac

            awk -v name="$name" -v threads="$threads" -v mode="$mode" '
                /^real / { real=$2 }
                /^user / { user=$2 }
                /^sys / { sys=$2 }
                END {
                    printf "%s\t%s\t%s\t%s\t%s\t%s\n",
                           name, threads, mode, real+0, user+0, sys+0
                }
            ' "$err" >> "$timings"
        done

        if cmp "$base_out" "$plan_out" >/dev/null 2>&1; then
            printf '%s\t%s\tbaseline_vs_plan\tok\n' "$name" "$threads" >> "$checks"
        else
            printf '%s\t%s\tbaseline_vs_plan\tDIFF\n' "$name" "$threads" >> "$checks"
        fi
        if [ "$keep_outputs" = 0 ]; then
            rm -f "$base_out" "$plan_out"
        fi
    done
done

printf 'wrote %s and %s\n' "$timings" "$checks"
