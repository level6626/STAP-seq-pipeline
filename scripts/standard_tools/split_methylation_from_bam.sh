#!/usr/bin/env bash
set -Eeuo pipefail

# Split an existing STAP aligned BAM by the R2 methylation code embedded in the
# read name, then make per-code read1 TSS BED/bedGraph/BigWig files.

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-STAP_TSS_500_oligos_S2}"
BAM="${BAM:-results/standard_tools/${SAMPLE}/${SAMPLE}.bowtie2.dedup.bam}"
OUTDIR="${OUTDIR:-results/standard_tools/${SAMPLE}/methylation_split}"
CHROM_SIZES="${CHROM_SIZES:-results/standard_tools/${SAMPLE}/reference/chrom.sizes}"
CANDIDATE_WINDOWS="${CANDIDATE_WINDOWS:-}"

# umi_tools extract appended R1_UMI + R2_17mer to the read name. The methylation
# code is the first 3 bases of R2, so it starts after the 8-bp R1 UMI.
APPENDED_R1_UMI_LENGTH="${APPENDED_R1_UMI_LENGTH:-8}"
INCLUDE_UNKNOWN="${INCLUDE_UNKNOWN:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${OUTDIR}" "${LOGDIR}"

trap 'echo "ERROR at line ${LINENO}. See ${LOGDIR} for logs." >&2' ERR

log_msg() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOGDIR}/methylation_split.progress.log" >&2
}

run_logged() {
    local name="$1"
    shift
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command:' -1
        printf ' %q' "$@"
        printf '\n'
        "$@"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

run_bash_logged() {
    local name="$1"
    local command="$2"
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command: %s\n' -1 "${command}"
        bash -o pipefail -c "${command}"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

require_file() {
    local path="$1"
    if [[ ! -s "${path}" ]]; then
        echo "Required file is missing or empty: ${path}" >&2
        exit 1
    fi
}

require_file "${BAM}"
require_file "${CHROM_SIZES}"

MAP_TSV="${OUTDIR}/methylation_code_outputs.tsv"
cat >"${MAP_TSV}" <<EOF
code	label	bam
TTT	100pct	${OUTDIR}/${SAMPLE}.TTT_100pct.bam
AAA	0pct	${OUTDIR}/${SAMPLE}.AAA_0pct.bam
CAT	60pct	${OUTDIR}/${SAMPLE}.CAT_60pct.bam
AGT	40pct	${OUTDIR}/${SAMPLE}.AGT_40pct.bam
TGA	20pct	${OUTDIR}/${SAMPLE}.TGA_20pct.bam
TAG	10pct	${OUTDIR}/${SAMPLE}.TAG_10pct.bam
CTA	1pct	${OUTDIR}/${SAMPLE}.CTA_1pct.bam
ATG	0p1pct	${OUTDIR}/${SAMPLE}.ATG_0p1pct.bam
EOF

if [[ "${INCLUDE_UNKNOWN}" == "1" || "${INCLUDE_UNKNOWN}" == "true" ]]; then
    printf 'unknown\tunknown\t%s/%s.unknown.bam\n' "${OUTDIR}" "${SAMPLE}" >>"${MAP_TSV}"
fi

SUMMARY_TSV="${OUTDIR}/${SAMPLE}.methylation_split.summary.tsv"
OFFSET=$((APPENDED_R1_UMI_LENGTH + 1))

run_bash_logged split_bam_by_methylation_code \
    "samtools view -h '${BAM}' | gawk -v map_tsv='${MAP_TSV}' -v summary_tsv='${SUMMARY_TSV}' -v threads='${THREADS}' -v offset='${OFFSET}' '
        BEGIN {
            FS = OFS = \"\t\"
            while ((getline line < map_tsv) > 0) {
                if (line !~ /^code\t/) {
                    split(line, f, \"\t\")
                    code = f[1]
                    label[code] = f[2]
                    bam[code] = f[3]
                    cmd[code] = \"samtools view -@ \" threads \" -b -o \\047\" f[3] \"\\047 -\"
                    order[++n] = code
                }
            }
            close(map_tsv)
        }
        /^@/ {
            for (i = 1; i <= n; i++) {
                print \$0 | cmd[order[i]]
            }
            next
        }
        {
            appended = \$1
            sub(/^.*_/, \"\", appended)
            code = \"unknown\"
            if (length(appended) >= offset + 2) {
                code = substr(appended, offset, 3)
            }
            if (!(code in cmd)) {
                code = \"unknown\"
            }
            if (code in cmd) {
                print \$0 | cmd[code]
                records[code]++
                if (and(\$2 + 0, 64)) {
                    read1[code]++
                }
            } else {
                skipped++
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                close(cmd[order[i]])
            }
            print \"code\", \"label\", \"bam\", \"bam_records\", \"read1_records\" > summary_tsv
            for (i = 1; i <= n; i++) {
                code = order[i]
                print code, label[code], bam[code], records[code] + 0, read1[code] + 0 >> summary_tsv
            }
            if (skipped > 0) {
                print \"skipped\", \"skipped\", \"NA\", skipped, \"NA\" >> summary_tsv
            }
        }
    '"

tail -n +2 "${MAP_TSV}" | while IFS=$'\t' read -r code label split_bam; do
    suffix="${code}_${label}"
    if [[ "${code}" == "${label}" ]]; then
        suffix="${code}"
    fi
    prefix="${OUTDIR}/${SAMPLE}.${suffix}"
    tss_bed="${prefix}.read1_tss.bed"
    tss_bedgraph="${prefix}.read1_tss.bedGraph"
    tss_bw="${prefix}.read1_tss.bw"

    run_logged "index_${suffix}" samtools index "${split_bam}"
    run_logged "flagstat_${suffix}" samtools flagstat "${split_bam}"

    run_bash_logged "tss_bed_${suffix}" \
        "samtools view -f 64 -F 2820 '${split_bam}' | gawk -f '${SCRIPT_DIR}/sam_read1_tss_to_bed.awk' | sort -k1,1 -k2,2n > '${tss_bed}'"

    run_bash_logged "tss_bedgraph_${suffix}" \
        "bedtools genomecov -bg -i '${tss_bed}' -g '${CHROM_SIZES}' | sort -k1,1 -k2,2n > '${tss_bedgraph}'"

    if [[ -s "${tss_bedgraph}" ]] && command -v bedGraphToBigWig >/dev/null 2>&1; then
        run_logged "tss_bigwig_${suffix}" bedGraphToBigWig "${tss_bedgraph}" "${CHROM_SIZES}" "${tss_bw}"
    fi

    if [[ -n "${CANDIDATE_WINDOWS}" ]]; then
        require_file "${CANDIDATE_WINDOWS}"
        run_bash_logged "candidate_counts_${suffix}" \
            "bedtools coverage -counts -a '${CANDIDATE_WINDOWS}' -b '${tss_bed}' > '${prefix}.candidate_window_tss_counts.tsv'"
    fi
done

log_msg "Methylation split complete for ${SAMPLE}"
column -t -s $'\t' "${SUMMARY_TSV}" | tee -a "${LOGDIR}/methylation_split.progress.log" >&2
