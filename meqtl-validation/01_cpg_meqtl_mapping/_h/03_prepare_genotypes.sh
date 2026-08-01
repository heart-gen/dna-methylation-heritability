#!/usr/bin/env bash
# Filter genotypes for a region-specific meQTL sample set.
# Requires plink2. Writes under meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/genotypes/
#
# Usage: 03_prepare_genotypes.sh <region> [AA|EA]
#   AA (default): TOPMed_LIBD.AA → meqtl_AA
#   EA: all_individuals TOPMed_LIBD → meqtl_EA
#
# Sample IDs: output FID=IID=BrNum so TensorQTL PgenReader (psam index = FID)
# matches phenotype/covariate BrNum columns.
set -euo pipefail

REGION="${1:?region required (caudate|dlpfc|hippocampus)}"
POPULATION="${2:-AA}"
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
MAP_ROOT="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping"
if [[ "${POPULATION}" == "AA" ]]; then
  PRE="${MAP_ROOT}/${REGION}/_m/preflight/sample_inclusion_primary.tsv"
  OUT_NAME="meqtl_AA"
  STAGE_NAME="TOPMed_LIBD.AA"
  SRC_PGEN="${ROOT}/inputs/genotypes/TOPMed_LIBD.AA.pgen"
  SRC_PVAR="${ROOT}/inputs/genotypes/TOPMed_LIBD.AA.pvar"
  SRC_PSAM="${ROOT}/inputs/genotypes/TOPMed_LIBD.AA.psam"
else
  PRE="${MAP_ROOT}/${REGION}/_m/preflight/${POPULATION}/sample_inclusion_primary.tsv"
  OUT_NAME="meqtl_${POPULATION}"
  STAGE_NAME="TOPMed_LIBD.all_individuals"
  SRC_PGEN="${ROOT}/inputs/genotypes/all_individuals/TOPMed_LIBD.pgen"
  SRC_PVAR="${ROOT}/inputs/genotypes/all_individuals/TOPMed_LIBD.pvar"
  SRC_PSAM="${ROOT}/inputs/genotypes/all_individuals/TOPMed_LIBD.psam"
fi
OUTDIR="${MAP_ROOT}/${REGION}/_m/genotypes"
STAGE="${MAP_ROOT}/_m/genotype_source"
mkdir -p "${OUTDIR}" "${STAGE}"

if [[ ! -f "${PRE}" ]]; then
  echo "Missing ${PRE}; run 00_preflight.py --population ${POPULATION} first" >&2
  exit 1
fi

if ! command -v plink2 >/dev/null 2>&1; then
  echo "plink2 not found on PATH" >&2
  exit 1
fi

GENO_PREFIX="${STAGE}/${STAGE_NAME}"

# Stage pfile with a plink2-compatible headered psam (do not modify source)
ln -sfn "${SRC_PGEN}" "${GENO_PREFIX}.pgen"
ln -sfn "${SRC_PVAR}" "${GENO_PREFIX}.pvar"
if [[ ! -f "${GENO_PREFIX}.psam" ]] || [[ "${SRC_PSAM}" -nt "${GENO_PREFIX}.psam" ]]; then
  {
    echo -e '#FID\tIID\tSEX'
    # Drop any existing header-like first line
    awk 'BEGIN{OFS="\t"} $1 !~ /^#/ {print $1,$2,(NF>=3?$3:"NA")}' "${SRC_PSAM}"
  } > "${GENO_PREFIX}.psam"
fi

# plink2 applies --update-ids before --keep, so KEEP must use post-update IDs
# (FID=IID=BrNum). IDMAP still uses original FID/IID from the staged psam.
KEEP="${OUTDIR}/keep_brnum_brnum.txt"
IDMAP="${OUTDIR}/update_ids_to_brnum.txt"
python3 - <<PY
from pathlib import Path

brnums = {
    line.strip().split("\t")[0]
    for line in Path("${PRE}").read_text().splitlines()[1:]
    if line.strip()
}
keep_rows = []
idmap_rows = []
matched = set()
for line in Path("${GENO_PREFIX}.psam").read_text().splitlines():
    parts = line.rstrip("\n").split("\t")
    if not parts or parts[0].startswith("#"):
        continue
    fid, iid = parts[0], parts[1] if len(parts) > 1 else parts[0]
    if fid in brnums:
        matched.add(fid)
        # after --update-ids, both FID and IID become BrNum
        keep_rows.append(f"{fid}\t{fid}")
        idmap_rows.append(f"{fid}\t{iid}\t{fid}\t{fid}")
missing = sorted(brnums - matched)
Path("${KEEP}").write_text("\n".join(keep_rows) + ("\n" if keep_rows else ""))
Path("${IDMAP}").write_text("\n".join(idmap_rows) + ("\n" if idmap_rows else ""))
print(f"keep n={len(keep_rows)}; missing_from_psam n={len(missing)}")
if missing:
    print("missing BrNums (first 20):", ",".join(missing[:20]))
if not keep_rows:
    raise SystemExit("No overlapping samples between inclusion list and psam")
PY

MAF=0.05
GENO_MISS=0.05
HWE=1e-6
TMP_PREFIX="${OUTDIR}/${OUT_NAME}.tmp"
OUT_PREFIX="${OUTDIR}/${OUT_NAME}"

# Note: TOPMed pvars may lack INFO/R2; imputation_r2_min from config cannot be applied.
plink2 \
  --pfile "${GENO_PREFIX}" \
  --update-ids "${IDMAP}" \
  --keep "${KEEP}" \
  --autosome \
  --maf "${MAF}" \
  --geno "${GENO_MISS}" \
  --hwe "${HWE}" \
  --make-pgen \
  --out "${TMP_PREFIX}"

for ext in pgen pvar psam; do
  mv -f "${TMP_PREFIX}.${ext}" "${OUT_PREFIX}.${ext}"
done
if [[ -f "${TMP_PREFIX}.log" ]]; then
  mv -f "${TMP_PREFIX}.log" "${OUT_PREFIX}.plink2.log"
fi

python3 - <<PY
from pathlib import Path
import datetime as dt

outdir = Path("${OUTDIR}")
out_name = "${OUT_NAME}"
psam_lines = (outdir / f"{out_name}.psam").read_text().splitlines()
n_samples = sum(1 for line in psam_lines if line.strip() and not line.startswith("#"))
# Confirm FID == IID == BrNum
bad = []
for line in psam_lines:
    if not line.strip() or line.startswith("#"):
        continue
    parts = line.split("\t")
    if parts[0] != parts[1]:
        bad.append(parts[0])
if bad:
    raise SystemExit(f"Expected FID=IID=BrNum; mismatches e.g. {bad[:5]}")
n_variants = sum(
    1 for line in (outdir / f"{out_name}.pvar").open()
    if line.strip() and not line.startswith("#")
)
summary = outdir / f"genotype_qc_summary_{out_name}.tsv"
header = "\t".join([
    "region", "population", "n_samples", "n_variants", "maf_min", "geno_missing_max",
    "hwe_p_min", "imputation_r2_min_applied", "imputation_r2_note",
    "sample_id_scheme", "source_pfile", "pfile_prefix", "generated_at",
])
row = "\t".join([
    "${REGION}",
    "${POPULATION}",
    str(n_samples),
    str(n_variants),
    "${MAF}",
    "${GENO_MISS}",
    "${HWE}",
    "NA",
    "INFO/R2 may be absent from source pvar; filter not applied",
    "FID=IID=BrNum",
    "${GENO_PREFIX}",
    str(outdir / out_name),
    dt.datetime.now().isoformat(timespec="seconds"),
])
text = header + "\n" + row + "\n"
summary.write_text(text)
if out_name == "meqtl_AA":
    (outdir / "genotype_qc_summary.tsv").write_text(text)
print(f"Wrote {summary}; samples={n_samples}; variants={n_variants}")
print(f"Wrote {outdir}/{out_name}.{{pgen,pvar,psam}}")
PY
