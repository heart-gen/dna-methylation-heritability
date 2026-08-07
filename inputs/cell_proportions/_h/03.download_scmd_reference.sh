#!/bin/bash
set -euo pipefail

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
REF_DIR="${ROOT}/inputs/cell_proportions/_m/reference"
SHA="9b4f52f4bc22ab1e39266f61a80233597e6b52c1"
EXPECTED="fd85dcaa013ccddbd4f564d03fdc3c64efcb0c59e4b6cf79156a1a626d1dc967"
URL="https://github.com/randel/scMD/archive/${SHA}.tar.gz"
ARCHIVE="${REF_DIR}/scMD-${SHA}.tar.gz"
PREFIX="scMD-${SHA}"
ANNO_VERSION="0.6.0"
ANNO_SHA="2c8128126b63e7fa805a5f3b02449367dca9c3be3eb5f6300acc718826590719"
ANNO_URL="https://bioconductor.org/packages/3.17/data/annotation/src/contrib/IlluminaHumanMethylationEPICanno.ilm10b4.hg19_${ANNO_VERSION}.tar.gz"
ANNO_ARCHIVE="${REF_DIR}/IlluminaHumanMethylationEPICanno.ilm10b4.hg19_${ANNO_VERSION}.tar.gz"

mkdir -p "${REF_DIR}"

if [[ ! -f "${ARCHIVE}" ]]; then
  curl -L --fail --retry 3 --max-time 600 -o "${ARCHIVE}.part" "${URL}"
  mv "${ARCHIVE}.part" "${ARCHIVE}"
fi

ACTUAL="$(sha256sum "${ARCHIVE}" | awk '{print $1}')"
if [[ "${ACTUAL}" != "${EXPECTED}" ]]; then
  echo "ERROR: scMD archive checksum mismatch: expected ${EXPECTED}, got ${ACTUAL}" >&2
  exit 1
fi

if [[ ! -f "${ANNO_ARCHIVE}" ]]; then
  curl -L --fail --retry 3 --max-time 600 -o "${ANNO_ARCHIVE}.part" "${ANNO_URL}"
  mv "${ANNO_ARCHIVE}.part" "${ANNO_ARCHIVE}"
fi
ANNO_ACTUAL="$(sha256sum "${ANNO_ARCHIVE}" | awk '{print $1}')"
if [[ "${ANNO_ACTUAL}" != "${ANNO_SHA}" ]]; then
  echo "ERROR: EPIC annotation checksum mismatch: expected ${ANNO_SHA}, got ${ANNO_ACTUAL}" >&2
  exit 1
fi

extract_one() {
  local member="$1"
  local output="$2"
  if [[ ! -f "${output}" ]]; then
    tar -xOzf "${ARCHIVE}" "${PREFIX}/${member}" > "${output}.part"
    mv "${output}.part" "${output}"
  fi
}

extract_one "data/Lee_7ct_WGBS.rda" "${REF_DIR}/Lee_7ct_WGBS.rda"
extract_one "data/Tian_7ct_WGBS.rda" "${REF_DIR}/Tian_7ct_WGBS.rda"
extract_one "data/Lee_7ct_450850.rda" "${REF_DIR}/Lee_7ct_450850.rda"
extract_one "data/Tian_7ct_450850.rda" "${REF_DIR}/Tian_7ct_450850.rda"
extract_one "Processed_data_450k850k/Lee_450k.csv" "${REF_DIR}/Lee_450k.csv"
extract_one "Processed_data_450k850k/Tian_450k.csv" "${REF_DIR}/Tian_450k.csv"
extract_one "Processed_data_450k850k/Lee_850k.csv" "${REF_DIR}/Lee_850k.csv"
extract_one "Processed_data_450k850k/Tian_850k.csv" "${REF_DIR}/Tian_850k.csv"

{
  printf 'asset\tsource_url\tversion\tcommit\tsha256\tgenome_build\n'
  printf 'scMD_source_archive\t%s\t1.0.0\t%s\t%s\tnot_applicable\n' "${URL}" "${SHA}" "${ACTUAL}"
  printf 'EPIC_hg19_annotation\t%s\t%s\tBioconductor_3.17\t%s\thg19\n' "${ANNO_URL}" "${ANNO_VERSION}" "${ANNO_ACTUAL}"
  printf 'Lee_7ct_WGBS\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Lee_7ct_WGBS.rda" | awk '{print $1}')"
  printf 'Tian_7ct_WGBS\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Tian_7ct_WGBS.rda" | awk '{print $1}')"
  printf 'Lee_7ct_450850\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Lee_7ct_450850.rda" | awk '{print $1}')"
  printf 'Tian_7ct_450850\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Tian_7ct_450850.rda" | awk '{print $1}')"
  printf 'Lee_450k_fallback\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Lee_450k.csv" | awk '{print $1}')"
  printf 'Tian_450k_fallback\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Tian_450k.csv" | awk '{print $1}')"
  printf 'Lee_850k_fallback\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Lee_850k.csv" | awk '{print $1}')"
  printf 'Tian_850k_fallback\t%s\t1.0.0\t%s\t%s\thg19\n' "${URL}" "${SHA}" "$(sha256sum "${REF_DIR}/Tian_850k.csv" | awk '{print $1}')"
} > "${REF_DIR}/scmd_reference_manifest.tsv"

echo "Prepared pinned scMD reference under ${REF_DIR}"
