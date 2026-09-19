#!/bin/bash
# 09 stage 17a -- genome-wide significant variants for every negative-control
# trait, plus PGC3 schizophrenia lifted to hg38.
#
#   ./17a_extract_gwas_leads.sh <out_dir>
#
# Streams each 7-9M-row file once with awk and keeps only p < 5e-8 rows, so
# the R stage never reads a whole GWAS. Output: <out_dir>/sig/{tag}.tsv with
# columns chrom pos pvalue variant_id, one file per trait, and
# <out_dir>/sig/PGC3_SCZ.tsv for the Module 09 trait under the same rule.
#
# Thresholds and paths come from config/gwas_negative_controls.yml, read here
# with a tiny python shim so this script and the R stage cannot disagree.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"
OUT=${1:?usage: 17a_extract_gwas_leads.sh <out_dir>}
CFG="${REPO_DIR}/config/gwas_negative_controls.yml"
mkdir -p "${OUT}/sig" "${OUT}/lift"

cfg () {  # dotted key -> scalar
    python3 - "$CFG" "$1" <<'EOF'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
for k in sys.argv[2].split("."):
    c = c[k]
print(c)
EOF
}
P=$(cfg locus_rule.genome_wide_p)
DIR=$(cfg collection.dir)
TPL=$(cfg collection.file_template)
C_CHR=$(cfg collection.columns.chrom); C_POS=$(cfg collection.columns.pos)
C_P=$(cfg collection.columns.pvalue);  C_ID=$(cfg collection.columns.variant_id)
META=$(cfg collection.metadata)

# Every file must carry the header the config's column indices assume.
EXPECT="variant_id	panel_variant_id	chromosome	position	effect_allele	non_effect_allele	current_build	frequency	sample_size	zscore	pvalue	effect_size	standard_error	imputation_status	n_cases"

n=0
for tag in $(cut -f2 "$META" | tail -n +2); do
    f="${DIR}/${TPL/\{tag\}/$tag}"
    if [ ! -f "$f" ]; then
        log_message "WARNING: no file for metadata tag ${tag}; skipped"
        continue
    fi
    # head closes the pipe early; under pipefail that SIGPIPE would abort the run.
    hdr=$( (zcat "$f" 2>/dev/null || true) | head -n 1)
    if [ "$hdr" != "$EXPECT" ]; then
        echo "ERROR: unexpected header in ${f}" >&2; exit 1
    fi
    out="${OUT}/sig/${tag}.tsv"
    printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${out}.tmp"
    zcat "$f" | awk -F'\t' -v p="$P" -v c="$C_CHR" -v s="$C_POS" -v q="$C_P" -v i="$C_ID" \
        'NR>1 && $q+0 < p+0 {print $c"\t"$s"\t"$q"\t"$i}' >> "${out}.tmp"
    mv "${out}.tmp" "$out"
    n=$((n+1))
    log_message "${tag}: $(($(wc -l < "$out") - 1)) genome-wide significant rows"
done
log_message "extracted ${n} traits"

# PGC3 SCZ: hg19 -> hg38 for the significant rows only.
SS=$(cfg scz_pgc3.sumstats_hg19)
CHAIN=$(cfg scz_pgc3.liftover_chain_hg19_to_hg38)
S_CHR=$(cfg scz_pgc3.columns.chrom); S_POS=$(cfg scz_pgc3.columns.pos)
S_P=$(cfg scz_pgc3.columns.pvalue);  S_ID=$(cfg scz_pgc3.columns.variant_id)
LIFT="${V2_LIFTOVER:-/projects/p32505/opt/envs/genomics/bin/liftOver}"
require_exec "$LIFT"
zcat "$SS" | grep -v '^##' | awk -F'\t' -v p="$P" -v c="$S_CHR" -v s="$S_POS" -v q="$S_P" -v i="$S_ID" \
    'NR>1 && $q+0 < p+0 {print "chr"$c"\t"$s-1"\t"$s"\t"$i"|"$q}' > "${OUT}/lift/pgc3.hg19.bed"
"$LIFT" "${OUT}/lift/pgc3.hg19.bed" "$CHAIN" "${OUT}/lift/pgc3.hg38.bed" "${OUT}/lift/pgc3.unmapped" >/dev/null 2>&1
printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${OUT}/sig/PGC3_SCZ.tsv.tmp"
awk -F'\t' '{split($4, a, "|"); print $1"\t"$3"\t"a[2]"\t"a[1]}' "${OUT}/lift/pgc3.hg38.bed" >> "${OUT}/sig/PGC3_SCZ.tsv.tmp"
mv "${OUT}/sig/PGC3_SCZ.tsv.tmp" "${OUT}/sig/PGC3_SCZ.tsv"
log_message "PGC3_SCZ: $(($(wc -l < "${OUT}/sig/PGC3_SCZ.tsv") - 1)) rows lifted, $(grep -c '^chr' "${OUT}/lift/pgc3.unmapped" || true) unmapped"
touch "${OUT}/sig/.complete"
