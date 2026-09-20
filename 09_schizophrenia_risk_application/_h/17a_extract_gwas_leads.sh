#!/bin/bash
# 09 stage 17a -- genome-wide significant variants for every negative-control
# trait, plus each focal GWAS (PGC3 schizophrenia and the multi-ancestry
# releases) lifted to hg38.
#
#   ./17a_extract_gwas_leads.sh <out_dir>
#
# Streams each 7-9M-row file once with awk and keeps only p < 5e-8 rows, so
# the R stage never reads a whole GWAS. Output: <out_dir>/sig/{tag}.tsv with
# columns chrom pos pvalue variant_id, one file per trait, including one per
# entry in the config's `focal_sumstats` list.
#
# Already-extracted tags are SKIPPED, so adding a focal trait costs one file
# rather than re-streaming the whole collection. FORCE=1 re-extracts everything.
#
# Thresholds and paths come from config/gwas_negative_controls.yml, read here
# with a tiny python shim so this script and the R stage cannot disagree.
#
# PARSING, and why it is not just `-F'\t'` (found 2026-09-19).
# PGC3_SCZ_wave3.primary carries a subset of rows whose first 12 values are
# SPACE-separated and then padded with empty tab fields. Under a plain tab split
# the p-value field of those rows is the empty string, and awk's `"" + 0 < 5e-8`
# is TRUE, so every malformed row entered the lead set as a spurious
# genome-wide hit. The parser below therefore
#   (1) requires the p-value to match a number before a row can pass, and
#   (2) re-splits a row on whitespace only when the tab split leaves a
#       non-numeric p-value, so well-formed rows are never re-interpreted;
#   (3) counts and reports rows it refused, rather than dropping them silently.
# A p-value of literal 0 is real underflow in several collection files and is
# kept. Verified after the fix: the PGC3 European extraction is unchanged.

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

focal_tsv () {  # one line per focal_sumstats entry
    python3 - "$CFG" <<'EOF'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
for e in c.get("focal_sumstats") or []:
    col = e["columns"]
    print("\t".join(str(x) for x in [
        e["tag"], e["path"], e.get("build", "hg38"),
        # NONE, never an empty field: bash collapses consecutive tabs even
        # when IFS is a tab alone, which would shift every later column.
        e.get("comment_prefix") or "NONE",
        col["chrom"], col["pos"], col["pvalue"], col["variant_id"]]))
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

# Shared row parser. -v mode=tsv emits chrom/pos/pvalue/variant_id;
# -v mode=bed emits a 4-column BED with "id|pvalue" in the name field.
read -r -d '' AWK_SIG <<'AWK' || true
function isnum(x) {
    return x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/
}
NR == 1 { next }
{
    cc = $c; ss = $s; qq = $q; ii = $i
    # Tab split left a non-numeric p-value: the row may be whitespace-separated
    # with tab padding. Re-split before judging it, then re-check.
    if (!isnum(qq)) {
        nf = split($0, f, /[ \t]+/)
        need = c; if (s > need) need = s; if (q > need) need = q; if (i > need) need = i
        if (nf >= need) { cc = f[c]; ss = f[s]; qq = f[q]; ii = f[i] }
    }
    if (!isnum(qq) || !isnum(ss)) { bad++; next }
    if (qq + 0 >= p + 0) next
    sub(/^chr/, "", cc)
    if (mode == "bed") print "chr" cc "\t" ss - 1 "\t" ss "\t" ii "|" qq
    else               print "chr" cc "\t" ss "\t" qq "\t" ii
}
END { if (bad) printf "MALFORMED_ROWS_SKIPPED=%d\n", bad > "/dev/stderr" }
AWK

# --------------------------------------------------------- collection traits
n=0; n_skip=0
for tag in $(cut -f2 "$META" | tail -n +2); do
    f="${DIR}/${TPL/\{tag\}/$tag}"
    if [ ! -f "$f" ]; then
        log_message "WARNING: no file for metadata tag ${tag}; skipped"
        continue
    fi
    if [ -s "${OUT}/sig/${tag}.tsv" ] && [ "${FORCE:-0}" != "1" ]; then
        n_skip=$((n_skip+1)); continue
    fi
    # head closes the pipe early; under pipefail that SIGPIPE would abort the run.
    hdr=$( (zcat "$f" 2>/dev/null || true) | head -n 1)
    if [ "$hdr" != "$EXPECT" ]; then
        echo "ERROR: unexpected header in ${f}" >&2; exit 1
    fi
    out="${OUT}/sig/${tag}.tsv"
    printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${out}.tmp"
    zcat "$f" | awk -F'\t' -v mode=tsv -v p="$P" -v c="$C_CHR" -v s="$C_POS" \
        -v q="$C_P" -v i="$C_ID" "$AWK_SIG" >> "${out}.tmp"
    mv "${out}.tmp" "$out"
    n=$((n+1))
    log_message "${tag}: $(($(wc -l < "$out") - 1)) genome-wide significant rows"
done
log_message "extracted ${n} collection traits (${n_skip} already present, skipped)"

# ------------------------------------------------------------------ focal GWAS
# Each entry of focal_sumstats: significant rows only, lifted to hg38 when the
# file is hg19. liftOver stderr is NOT discarded -- it is the only thing that
# reports a malformed BED, and discarding it hid exactly that once already.
CHAIN=$(cfg liftover_chain_hg19_to_hg38)
LIFT="${V2_LIFTOVER:-/projects/p32505/opt/envs/genomics/bin/liftOver}"
require_exec "$LIFT"

while IFS=$'\t' read -r tag path build comment c_chr c_pos c_p c_id; do
    [ -n "${tag:-}" ] || continue
    if [ ! -f "$path" ]; then
        echo "ERROR: focal sumstats missing: ${path}" >&2; exit 1
    fi
    if [ -s "${OUT}/sig/${tag}.tsv" ] && [ "${FORCE:-0}" != "1" ]; then
        log_message "${tag}: already present, skipped"; continue
    fi
    [ "$comment" = "NONE" ] && comment=""
    strip () { if [ -n "$comment" ]; then grep -v "^${comment}"; else cat; fi; }

    if [ "$build" = "hg38" ]; then
        printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${OUT}/sig/${tag}.tsv.tmp"
        (zcat "$path" 2>/dev/null || true) | strip \
            | awk -F'\t' -v mode=tsv -v p="$P" -v c="$c_chr" -v s="$c_pos" \
                -v q="$c_p" -v i="$c_id" "$AWK_SIG" >> "${OUT}/sig/${tag}.tsv.tmp"
        mv "${OUT}/sig/${tag}.tsv.tmp" "${OUT}/sig/${tag}.tsv"
        log_message "${tag}: $(($(wc -l < "${OUT}/sig/${tag}.tsv") - 1)) rows (already hg38)"
        continue
    fi

    bed19="${OUT}/lift/${tag}.hg19.bed"
    (zcat "$path" 2>/dev/null || true) | strip \
        | awk -F'\t' -v mode=bed -v p="$P" -v c="$c_chr" -v s="$c_pos" \
            -v q="$c_p" -v i="$c_id" "$AWK_SIG" > "$bed19"
    if [ ! -s "$bed19" ]; then
        # A real and reportable outcome: an underpowered ancestry-specific
        # release can carry no genome-wide significant variant at all.
        printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${OUT}/sig/${tag}.tsv"
        log_message "${tag}: 0 genome-wide significant variants; empty lead set written"
        continue
    fi
    "$LIFT" "$bed19" "$CHAIN" "${OUT}/lift/${tag}.hg38.bed" "${OUT}/lift/${tag}.unmapped" \
        > "${OUT}/lift/${tag}.liftover.log" 2>&1 \
        || { echo "ERROR: liftOver failed for ${tag}; see ${OUT}/lift/${tag}.liftover.log" >&2
             tail -3 "${OUT}/lift/${tag}.liftover.log" >&2; exit 1; }
    printf 'chrom\tpos\tpvalue\tvariant_id\n' > "${OUT}/sig/${tag}.tsv.tmp"
    awk -F'\t' '{split($4, a, "|"); print $1"\t"$3"\t"a[2]"\t"a[1]}' "${OUT}/lift/${tag}.hg38.bed" \
        >> "${OUT}/sig/${tag}.tsv.tmp"
    mv "${OUT}/sig/${tag}.tsv.tmp" "${OUT}/sig/${tag}.tsv"
    log_message "${tag}: $(($(wc -l < "${OUT}/sig/${tag}.tsv") - 1)) rows lifted, $(grep -c '^chr' "${OUT}/lift/${tag}.unmapped" || true) unmapped"
done < <(focal_tsv)

touch "${OUT}/sig/.complete"
