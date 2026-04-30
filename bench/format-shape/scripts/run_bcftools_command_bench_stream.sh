#!/bin/bash
set -euo pipefail

# Streaming variant of run_bcftools_command_bench.sh for very large VCFs.
# It runs the same command families but pipes command output through cksum,
# avoiding temporary BCF/text outputs that can be hundreds of GiB on full
# cohort chromosomes.  Baseline and planned checksums are compared.

bcftools=${BCFTOOLS:-bcftools}
inputs=${1:-bench/format-shape/large/bcftools-command-inputs.tsv}
outdir=${OUTDIR:-bench/format-shape/large/results-bcftools-commands-stream}
threads_list=${THREADS_LIST:-0}
commands=${COMMANDS:-view_bcf view_sites query_sites query_format stats filter_gt}
query_sample_count=${QUERY_SAMPLE_COUNT:-2}
mkdir -p "$outdir"

timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"
cmds_out="$outdir/commands.tsv"
checksums="$outdir/checksums.tsv"

printf 'name\tcommand\tthreads\tmode\treal\tuser\tsys\n' > "$timings"
printf 'name\tcommand\tthreads\tcomparison\tstatus\n' > "$checks"
printf 'name\tcommand\tthreads\tmode\tcksum\tbytes\n' > "$checksums"
printf 'command\tdescription\n' > "$cmds_out"
printf 'view_bcf\tbcftools view --no-version -Ob -l 0 streamed to cksum\n' >> "$cmds_out"
printf 'view_sites\tbcftools view --no-version -G -Ob -l 0 streamed to cksum\n' >> "$cmds_out"
printf 'query_sites\tbcftools query fixed site fields streamed to cksum\n' >> "$cmds_out"
printf 'query_format\tbcftools query GT for first QUERY_SAMPLE_COUNT samples streamed to cksum\n' >> "$cmds_out"
printf 'stats\tbcftools stats streamed to cksum\n' >> "$cmds_out"
printf 'filter_gt\tbcftools view -i GT="alt" for first QUERY_SAMPLE_COUNT samples streamed to cksum\n' >> "$cmds_out"
printf 'merge_self\tbcftools merge --no-index --force-samples streamed to cksum\n' >> "$cmds_out"

run_one()
{
    local mode=$1
    local command=$2
    local threads=$3
    local path=$4
    local sum_out=$5
    local err=$6
    local sample_args=$7
    local plan=0
    local thread_args=

    if [ "$mode" = plan ]; then
        plan=1
    fi
    if [ "$threads" != 0 ]; then
        thread_args="--threads $threads"
    fi

    case "$command" in
        view_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 $thread_args \
                -o - "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        view_sites)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -G -Ob -l 0 $thread_args \
                -o - "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        query_sites)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
                "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        query_format)
            # shellcheck disable=SC2086
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query $sample_args -f '%CHROM\t%POS[\t%GT]\n' \
                "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        stats)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" stats "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        filter_gt)
            # shellcheck disable=SC2086
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 $thread_args \
                $sample_args -i 'GT="alt"' -o - "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        merge_self)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" merge --no-index --force-samples --no-version -Ob \
                $thread_args -o - "$path" "$path" 2> "$err" | cksum > "$sum_out"
            ;;
        *)
            printf 'unknown command: %s\n' "$command" >&2
            return 1
            ;;
    esac
}

tail -n +2 "$inputs" | while IFS=$'\t' read -r name path source
do
    samples=$("$bcftools" query -l "$path" | awk -v n="$query_sample_count" '
        NR <= n { if (s) s = s "," $0; else s = $0 }
        END { print s }
    ')
    sample_args=
    if [ -n "$samples" ]; then
        sample_args="-s $samples"
    fi

    for command in $commands
    do
        case "$command" in
            query_format|filter_gt|merge_self)
                if [ -z "$sample_args" ]; then
                    for threads in $threads_list
                    do
                        printf '%s\t%s\t%s\tbaseline_vs_plan\tskipped_no_samples\n' \
                               "$name" "$command" "$threads" >> "$checks"
                    done
                    continue
                fi
                ;;
        esac

        for threads in $threads_list
        do
            base_sum="$outdir/$name.$command.t$threads.baseline.cksum"
            plan_sum="$outdir/$name.$command.t$threads.plan.cksum"

            for mode in baseline plan
            do
                err="$outdir/$name.$command.t$threads.$mode.stderr"
                sum="$outdir/$name.$command.t$threads.$mode.cksum"
                run_one "$mode" "$command" "$threads" "$path" "$sum" "$err" "$sample_args"

                awk -v name="$name" -v command="$command" \
                    -v threads="$threads" -v mode="$mode" '
                    /^real / { real=$2 }
                    /^user / { user=$2 }
                    /^sys / { sys=$2 }
                    END {
                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                               name, command, threads, mode,
                               real+0, user+0, sys+0
                    }
                ' "$err" >> "$timings"

                awk -v name="$name" -v command="$command" \
                    -v threads="$threads" -v mode="$mode" '
                    { printf "%s\t%s\t%s\t%s\t%s\t%s\n",
                             name, command, threads, mode, $1, $2 }
                ' "$sum" >> "$checksums"
            done

            if cmp "$base_sum" "$plan_sum" >/dev/null 2>&1; then
                printf '%s\t%s\t%s\tbaseline_vs_plan\tok\n' \
                       "$name" "$command" "$threads" >> "$checks"
            else
                printf '%s\t%s\t%s\tbaseline_vs_plan\tDIFF\n' \
                       "$name" "$command" "$threads" >> "$checks"
            fi
        done
    done
done

printf 'wrote %s, %s, %s, and %s\n' "$timings" "$checks" "$checksums" "$cmds_out"
