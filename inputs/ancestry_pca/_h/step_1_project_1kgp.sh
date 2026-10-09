#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=ancestry_pca
#SBATCH --output=logs/ancestry_pca.%A.log
#
# inputs/ancestry_pca step 1: place the cohort in the 1000 Genomes PC space.
#
# The ancestry-confirmation supplement used to overlay the cohort's own
# genotype PCs on a separately computed 1000 Genomes PCA. Two PCAs have
# unrelated axes, so the overlay put white American donors beside AFR. Here
# one PCA is fit on 1000 Genomes alone, and both the reference panel and the
# cohort are scored with the same allele weights, so every point shares a basis.
#
# Variants: autosomal, biallelic, A/C/G/T SNPs with MAF >= 0.05 and missingness
# <= 0.05 in both panels, keyed on chr:pos, with identical allele sets and no
# strand-ambiguous (A/T, C/G) pairs, then LD-pruned in the reference.
#
# Usage, from inputs/ancestry_pca/_m:
#   mkdir -p logs && sbatch ../_h/step_1_project_1kgp.sh
#
# Output (small, tracked): projected_pcs.tsv (source, group-free id, PC1-PC10),
# variant_counts.tsv, 1kgp_shared.eigenval. Donor IDs are BrNum; the array
# barcode is dropped.

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

COHORT="$REPO_DIR/inputs/genotypes/all_individuals/TOPMed_LIBD"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/_m/1kGP"
OUT="$REPO_DIR/inputs/ancestry_pca/_m"
TMP="$OUT/tmp"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

require_exec "$PLINK2"
for f in "$COHORT.pgen" "$REF.pgen"; do require_file "$f"; done
log_job_info
mkdir -p "$TMP"

# Shared filter, applied identically to both panels. --set-all-var-ids is
# applied at load, so --extract below matches the chr:pos IDs.
FILT=(--autosome --snps-only just-acgt --max-alleles 2 --maf 0.05 --geno 0.05
      --set-all-var-ids '@:#' --rm-dup exclude-all --threads "$THREADS")

log_message "**** Candidate variants in each panel ****"
"$PLINK2" --pfile "$COHORT" "${FILT[@]}" --make-just-pvar --out "$TMP/cohort_all"
"$PLINK2" --pfile "$REF"    "${FILT[@]}" --make-just-pvar --out "$TMP/ref_all"

log_message "**** Shared chr:pos with identical, strand-unambiguous allele sets ****"
awk 'BEGIN{FS=OFS="\t"}
     function key(a,b){ return (a<b) ? a"/"b : b"/"a }
     function ambig(a,b){ return (a b=="AT"||a b=="TA"||a b=="CG"||a b=="GC") }
     /^#/ {next}
     FNR==NR { c[$3]=key($4,$5); next }
     ($3 in c) && c[$3]==key($4,$5) && !ambig($4,$5) { print $3 }' \
    "$TMP/cohort_all.pvar" "$TMP/ref_all.pvar" > "$TMP/shared.ids"

log_message "**** LD-prune in the reference ****"
"$PLINK2" --pfile "$REF" "${FILT[@]}" --extract "$TMP/shared.ids" \
    --indep-pairwise 500 50 0.1 --out "$TMP/ref_prune"

log_message "**** Reference PCA with allele weights ****"
"$PLINK2" --pfile "$REF" "${FILT[@]}" --extract "$TMP/ref_prune.prune.in" \
    --freq counts --pca 10 allele-wts --out "$TMP/1kgp_shared"

# Locate the A1 and PC1..PC10 columns by name; their positions vary by build.
HDR=$(head -n 1 "$TMP/1kgp_shared.eigenvec.allele")
A1_COL=$(echo "$HDR" | tr '\t' '\n' | grep -n -x 'A1' | cut -d: -f1)
PC1_COL=$(echo "$HDR" | tr '\t' '\n' | grep -n -x 'PC1' | cut -d: -f1)
[ -n "$A1_COL" ] && [ -n "$PC1_COL" ] || { echo "ERROR: unexpected eigenvec.allele header: $HDR" >&2; exit 1; }

log_message "**** Score both panels with the same weights ****"
for PANEL in ref cohort; do
    SRC="$REF"; [ "$PANEL" = "cohort" ] && SRC="$COHORT"
    "$PLINK2" --pfile "$SRC" "${FILT[@]}" --extract "$TMP/ref_prune.prune.in" \
        --read-freq "$TMP/1kgp_shared.acount" \
        --score "$TMP/1kgp_shared.eigenvec.allele" 2 "$A1_COL" header-read \
                no-mean-imputation variance-standardize \
        --score-col-nums "${PC1_COL}-$((PC1_COL + 9))" \
        --out "$TMP/${PANEL}_projected"
done

log_message "**** Write tracked outputs ****"
# The cohort .sscore carries FID = BrNum and IID = array barcode; keep BrNum.
# The reference has no FID column, so IID is the 1000 Genomes sample ID.
emit() {
    awk -v src="$2" 'BEGIN{FS=OFS="\t"}
         NR==1 { for (i=1;i<=NF;i++) { if ($i=="#FID") f=i; if ($i=="IID"||$i=="#IID") d=i;
                                        if ($i ~ /^PC[0-9]+_AVG$/) pc[++n]=i }
                 next }
         { id = f ? $f : $d; line = src OFS id
           for (k=1;k<=n;k++) line = line OFS $pc[k]
           print line }' "$1"
}
{
    printf 'source\tid'; for k in $(seq 1 10); do printf '\tPC%s' "$k"; done; printf '\n'
    emit "$TMP/ref_projected.sscore" "1000 Genomes"
    emit "$TMP/cohort_projected.sscore" "this study"
} > "$OUT/projected_pcs.tsv.tmp"
mv "$OUT/projected_pcs.tsv.tmp" "$OUT/projected_pcs.tsv"
cp "$TMP/1kgp_shared.eigenval" "$OUT/1kgp_shared.eigenval"

{
    printf 'step\tn_variants\n'
    printf 'cohort_candidates\t%s\n' "$(grep -vc '^#' "$TMP/cohort_all.pvar")"
    printf 'reference_candidates\t%s\n' "$(grep -vc '^#' "$TMP/ref_all.pvar")"
    printf 'shared_allele_matched\t%s\n' "$(wc -l < "$TMP/shared.ids")"
    printf 'ld_pruned_scored\t%s\n' "$(wc -l < "$TMP/ref_prune.prune.in")"
} > "$OUT/variant_counts.tsv"

sha256sum "$OUT/projected_pcs.tsv" "$OUT/1kgp_shared.eigenval" "$OUT/variant_counts.tsv" \
    > "$OUT/output_checksums.tsv"
rm -rf "$TMP"
log_message "**** Job ends ****"
