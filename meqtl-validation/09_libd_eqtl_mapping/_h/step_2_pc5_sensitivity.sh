#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=libd_pc5_tqtl
#SBATCH --output=logs/libd_pc5_tqtl.%j.log

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
MAPPER="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h/04_tensorqtl_map.py"
REGION="${REGION:-caudate}"
OUT="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/${REGION}/genes_rpkm_pc5"
PREP="${OUT}/prepared"
STD="${OUT}/standard"
TQTL="${OUT}/tensorqtl"
PHENO="${PREP}/genes.expression.bed.gz"
COV="${STD}/covariates.txt"
GENO="${ROOT}/inputs/genotypes/TOPMed_LIBD.AA"
PREFIX="libd_aa_${REGION}_genes_rpkm_pc5"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

mkdir -p "${TQTL}"
head -1 "${COV}"
python "${MAPPER}" \
  --region "${REGION}" \
  --mode cis \
  --phenotype-bed "${PHENO}" \
  --covariates "${COV}" \
  --genotype-prefix "${GENO}" \
  --outdir "${TQTL}" \
  --prefix "${PREFIX}" \
  --window 500000 \
  --maf 0.01 \
  --device "${DEVICE:-gpu}" \
  --chunk-size "${CHUNK_SIZE:-chr}" \
  --seed 13131313 \
  --fdr 0.05

# Append to comparison table
python - <<'PY'
from pathlib import Path
import pandas as pd
ROOT=Path("/projects/b1213/users/kynon/projects/dna-methylation-heritability")
M=ROOT/"meqtl-validation/09_libd_eqtl_mapping/_m/caudate"
cis=M/"genes_rpkm_pc5/tensorqtl/libd_aa_caudate_genes_rpkm_pc5.cis_qtl.txt.gz"
cov=M/"genes_rpkm_pc5/standard/covariates.txt"
prep=M/"genes_rpkm/prepared/prep_summary.tsv"
df=pd.read_csv(cis,sep="\t")
cov_cols=pd.read_csv(cov,sep="\t",nrows=0).columns
n_pcs=sum(c.startswith("PC") and c[2:].isdigit() for c in cov_cols)
prep_df=pd.read_csv(prep,sep="\t")
row={
 "label":"rpkm_mean0.2_pc5",
 "cis_path":str(cis),
 "n_phenotypes_mapped":len(df),
 "n_egene_fdr05":int((df.qval<=0.05).sum()),
 "n_egene_fdr10":int((df.qval<=0.10).sum()),
 "n_egene_fdr20":int((df.qval<=0.20).sum()),
 "min_qval":float(df.qval.min()),
 "n_expr_pcs":n_pcs,
 "n_samples":int(prep_df.iloc[0].n_samples),
 "n_features":int(prep_df.iloc[0].n_features),
 "phenotype":prep_df.iloc[0].phenotype,
 "feature_filter":prep_df.iloc[0].feature_filter,
 "norm_method":prep_df.iloc[0].norm_method,
 "status":"ok",
}
out=M/"sensitivity_rpkm_vs_cpm.tsv"
base=pd.read_csv(out,sep="\t") if out.exists() else pd.DataFrame()
base=base[base.label!="rpkm_mean0.2_pc5"]
base=pd.concat([base,pd.DataFrame([row])],ignore_index=True)
base.to_csv(out,sep="\t",index=False)
print(base.to_string(index=False))
PY

log_message "**** Job ends ****"
