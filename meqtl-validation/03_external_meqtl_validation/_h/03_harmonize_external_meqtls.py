#!/usr/bin/env python3
"""Harmonize public brain cis-meQTL catalogs to hg38 and annotate VMR CpGs.

Produces:
  _m/harmonized/{resource_id}.hg38.tsv.gz
  _m/harmonized/{resource_id}.overlap_summary.tsv
  _m/harmonized/{resource_id}.{region}.vmr_support.tsv.gz

Columns for the primary harmonized table (CpG-level):
  resource_id, tissue_region, phenotype_id, chrom, pos_1based, probe_id,
  genome_build_source, genome_build_harmonized, lift_status,
  external_pvalue, external_fdr, external_beta, lead_snp, snp_chrom, snp_pos,
  external_meqtl_support

VMR-annotated per-region tables add vmr_id for Phase 3 modeling.

BrainSeq full genome-wide catalogs require Synapse syn25992404. Until that
token is valid, this script harmonizes:
  - Nature SCZ-risk SNP–CpG tables as brainseq_wgbs_meqtl_scz_subset
  - GEO Jaffe DLPFC 450K all-pairs (jaffe_dlpfc_450k_meqtl)
  - Schulz hippocampus array meQTLs (schulz_hippocampus_array_meqtl)
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _lib.io_utils import PROJECT_ROOT, write_tsv  # noqa: E402

BASE = PROJECT_ROOT / "meqtl-validation" / "03_external_meqtl_validation" / "_m"
RAW = BASE / "raw"
HARM = BASE / "harmonized"
SUPPORT = BASE / "support"
PHASE1 = PROJECT_ROOT / "meqtl-validation" / "01_cpg_meqtl_mapping"
REGIONS = ("caudate", "dlpfc", "hippocampus")
CHAIN = SUPPORT / "hg19ToHg38.over.chain.gz"
LIFTOVER = Path("/projects/p32505/opt/envs/genomics/bin/liftOver")
DEFAULT_FDR = 0.05

# Illumina 450K assayed-CpG universe. Both Jaffe and Schulz are 450K studies, so the
# array manifest -- not the published results table -- defines which CpGs were
# interrogated. A manifest probe absent from a results table is a tested negative.
ANNO_RSCRIPT = Path("/projects/p32505/opt/envs/epigenomics/bin/Rscript")
ANNO_PACKAGE = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
UNIVERSE_450K = SUPPORT / "450k_universe_hg38.tsv.gz"
# Guard against silently joining a results table to the wrong universe.
MIN_SIG_PROBE_COVERAGE = 0.95


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--resources",
        nargs="+",
        default=["jaffe_dlpfc_450k_meqtl", "schulz_hippocampus_array_meqtl", "brainseq_wgbs_meqtl_scz_subset"],
    )
    p.add_argument("--fdr", type=float, default=DEFAULT_FDR)
    p.add_argument("--regions", nargs="+", default=list(REGIONS))
    return p.parse_args()


def ensure_chain() -> Path:
    if CHAIN.exists() and CHAIN.stat().st_size > 0:
        return CHAIN
    SUPPORT.mkdir(parents=True, exist_ok=True)
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    subprocess.check_call(["curl", "-L", "--fail", "-o", str(CHAIN), url])
    return CHAIN


def lift_positions(df: pd.DataFrame, chrom_col: str, pos_col: str) -> pd.DataFrame:
    """Lift unique 1-based CpG positions from hg19 to hg38."""
    ensure_chain()
    uniq = (
        df[[chrom_col, pos_col]]
        .dropna()
        .astype({pos_col: int})
        .drop_duplicates()
        .copy()
    )
    uniq["chrom"] = uniq[chrom_col].astype(str).map(lambda c: c if str(c).startswith("chr") else f"chr{c}")
    uniq["start0"] = uniq[pos_col].astype(int) - 1
    uniq["end0"] = uniq[pos_col].astype(int)
    uniq["uid"] = [f"u{i}" for i in range(len(uniq))]

    with tempfile.TemporaryDirectory(prefix="meqtl_lift_") as tmp:
        tmp_path = Path(tmp)
        bed_in = tmp_path / "in.bed"
        bed_out = tmp_path / "out.bed"
        bed_unmapped = tmp_path / "unmapped.bed"
        with bed_in.open("w") as handle:
            for row in uniq.itertuples(index=False):
                handle.write(f"{row.chrom}\t{row.start0}\t{row.end0}\t{row.uid}\n")
        cmd = [str(LIFTOVER), str(bed_in), str(CHAIN), str(bed_out), str(bed_unmapped)]
        subprocess.check_call(cmd)

        mapped = {}
        if bed_out.exists():
            with bed_out.open() as handle:
                for line in handle:
                    chrom, start, end, uid = line.rstrip("\n").split("\t")[:4]
                    mapped[uid] = (chrom, int(end))  # 1-based end for CpG
        unmapped = set()
        if bed_unmapped.exists():
            with bed_unmapped.open() as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 4:
                        unmapped.add(parts[3])

    uniq["lift_status"] = uniq["uid"].map(
        lambda u: "mapped" if u in mapped else ("unmapped" if u in unmapped else "missing")
    )
    uniq["chrom_hg38"] = uniq["uid"].map(lambda u: mapped[u][0] if u in mapped else pd.NA)
    uniq["pos_hg38"] = uniq["uid"].map(lambda u: mapped[u][1] if u in mapped else pd.NA)
    return uniq[[chrom_col, pos_col, "chrom_hg38", "pos_hg38", "lift_status"]]


def ensure_450k_universe() -> pd.DataFrame:
    """Return the Illumina 450K assayed-CpG universe lifted to hg38.

    Cached at `_m/support/450k_universe_hg38.tsv.gz`. Columns: probe_id, chrom,
    pos_1based. Only probes that lift cleanly are retained; an unlifted probe cannot
    be matched to a WGBS VMR CpG and so cannot serve as a negative either.
    """
    if UNIVERSE_450K.exists() and UNIVERSE_450K.stat().st_size > 0:
        u = pd.read_csv(UNIVERSE_450K, sep="\t", dtype={"probe_id": str, "chrom": str})
        u["pos_1based"] = u["pos_1based"].astype(int)
        print(f"450K universe (cached): {len(u):,} probes")
        return u

    SUPPORT.mkdir(parents=True, exist_ok=True)
    if not ANNO_RSCRIPT.exists():
        raise SystemExit(f"Rscript not found for manifest extraction: {ANNO_RSCRIPT}")
    print(f"Building 450K assayed universe from {ANNO_PACKAGE}")

    with tempfile.TemporaryDirectory(prefix="meqtl_450k_") as tmp:
        tmp_path = Path(tmp)
        bed_hg19 = tmp_path / "450k_hg19.bed"
        bed_hg38 = tmp_path / "450k_hg38.bed"
        bed_unmapped = tmp_path / "450k_unmapped.bed"
        r_code = f"""
suppressMessages(library({ANNO_PACKAGE}))
loc <- getAnnotation({ANNO_PACKAGE})
df <- data.frame(
  chr   = as.character(loc$chr),
  start = as.integer(loc$pos) - 1L,
  end   = as.integer(loc$pos),
  name  = rownames(loc),
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$start) & !is.na(df$chr), ]
write.table(df, "{bed_hg19}", sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("manifest probes:", nrow(df), "\\n")
"""
        subprocess.check_call([str(ANNO_RSCRIPT), "-e", r_code])
        ensure_chain()
        subprocess.check_call(
            [str(LIFTOVER), str(bed_hg19), str(CHAIN), str(bed_hg38), str(bed_unmapped)]
        )
        n_source = sum(1 for _ in bed_hg19.open())
        u = pd.read_csv(
            bed_hg38, sep="\t", header=None,
            names=["chrom", "start0", "pos_1based", "probe_id"],
            dtype={"chrom": str, "probe_id": str},
        )

    u = u[["probe_id", "chrom", "pos_1based"]].drop_duplicates(subset=["probe_id"], keep="first")
    u["pos_1based"] = u["pos_1based"].astype(int)
    u.to_csv(UNIVERSE_450K, sep="\t", index=False, compression="gzip")
    print(
        f"450K universe: {len(u):,}/{n_source:,} probes lifted to hg38 "
        f"({len(u) / n_source:.4%}); cached at {UNIVERSE_450K}"
    )
    return u


def build_array_resource(
    resource_id: str,
    tissue_region: str,
    sig: pd.DataFrame,
    fdr: float,
) -> tuple[pd.DataFrame, dict]:
    """Expand an array meQTL results table onto the full 450K assayed universe.

    `sig` must be one row per significant probe with columns
    probe_id, external_pvalue, external_fdr, external_beta, lead_snp, snp_chrom, snp_pos.
    Every manifest probe becomes a row; probes absent from `sig` are tested negatives.
    """
    universe = ensure_450k_universe()
    sig = sig.drop_duplicates(subset=["probe_id"], keep="first").copy()

    n_sig = len(sig)
    n_matched = int(sig["probe_id"].isin(set(universe["probe_id"])).sum())
    coverage = n_matched / n_sig if n_sig else 0.0
    print(
        f"  {resource_id}: {n_matched:,}/{n_sig:,} significant probes on the 450K "
        f"manifest ({coverage:.1%})"
    )
    if coverage < MIN_SIG_PROBE_COVERAGE:
        raise SystemExit(
            f"{resource_id}: only {coverage:.1%} of significant probes map to the 450K "
            f"manifest (threshold {MIN_SIG_PROBE_COVERAGE:.0%}). This resource is "
            "probably not a 450K study; do not manufacture negatives for it."
        )

    out = universe.merge(sig, on="probe_id", how="left")
    supported = out["external_fdr"].notna() & (out["external_fdr"] <= fdr)
    out["resource_id"] = resource_id
    out["tissue_region"] = tissue_region
    out["phenotype_id"] = out["chrom"].astype(str) + "_" + out["pos_1based"].astype(int).astype(str)
    out["genome_build_source"] = "hg19"
    out["genome_build_harmonized"] = "hg38"
    out["lift_status"] = "mapped"
    out["external_assayed"] = 1
    out["external_meqtl_support"] = supported.astype(int)
    out["assay_universe_complete"] = True
    out["assay_universe_source"] = f"{ANNO_PACKAGE}_manifest"

    stats = {
        "n_probes_universe": int(len(out)),
        "n_probes_supported": int(supported.sum()),
        "n_sig_probes_source": n_sig,
        "n_sig_probes_on_manifest": n_matched,
        "sig_probe_manifest_coverage": round(coverage, 6),
        "assay_universe_source": f"{ANNO_PACKAGE}_manifest",
        "fdr_threshold": fdr,
    }
    print(
        f"  {resource_id}: universe {stats['n_probes_universe']:,} probes; "
        f"{stats['n_probes_supported']:,} supported "
        f"({stats['n_probes_supported'] / max(stats['n_probes_universe'], 1):.3%})"
    )
    return out, stats


def load_cpg_vmr_maps(regions: list[str]) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for region in regions:
        prep = PHASE1 / region / "_m" / "prepared"
        files = sorted(prep.glob("cpg_vmr_map.chr*.tsv"))
        if not files:
            print(f"WARNING: no cpg_vmr_map for {region}")
            continue
        frames = [pd.read_csv(f, sep="\t", dtype={"phenotype_id": str, "vmr_id": str}) for f in files]
        m = pd.concat(frames, ignore_index=True)
        m["chrom"] = m["chrom"].astype(str)
        m["pos_1based"] = m["pos_1based"].astype(int)
        out[region] = m[["phenotype_id", "chrom", "pos_1based", "vmr_id"]].drop_duplicates()
        print(f"{region}: {len(out[region]):,} VMR CpGs")
    return out


def annotate_regions(
    cpg_level: pd.DataFrame,
    resource_id: str,
    region_maps: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Restrict to CpGs demonstrably assayed by the external resource.

    A WGBS CpG absent from an array catalog is unobserved, not a negative. The
    harmonized input must therefore carry the complete assayed-CpG universe.
    """
    if "external_assayed" not in cpg_level.columns:
        raise SystemExit(f"{resource_id}: harmonized table lacks external_assayed")
    assay = cpg_level.loc[cpg_level["external_assayed"] == 1].copy()
    assay = assay.dropna(subset=["chrom", "pos_1based"])
    assay["chrom"] = assay["chrom"].astype(str)
    assay["pos_1based"] = assay["pos_1based"].astype(int)
    if "external_fdr" in assay.columns:
        assay = assay.sort_values(
            ["external_meqtl_support", "external_fdr"], ascending=[False, True]
        )
    assay = assay.drop_duplicates(subset=["chrom", "pos_1based"], keep="first")

    summaries = []
    for region, mmap in region_maps.items():
        merged = mmap.merge(
            assay[["chrom", "pos_1based", "external_assayed", "external_meqtl_support", "external_fdr", "external_pvalue", "lead_snp", "probe_id", "assay_universe_complete"]],
            on=["chrom", "pos_1based"],
            how="inner",
        )
        merged["external_meqtl_support"] = merged["external_meqtl_support"].astype(int)
        merged["resource_id"] = resource_id
        merged["region"] = region
        out = HARM / f"{resource_id}.{region}.vmr_support.tsv.gz"
        merged.to_csv(out, sep="\t", index=False, compression="gzip")
        n_support = int(merged["external_meqtl_support"].sum())
        summaries.append({
            "resource_id": resource_id,
            "region": region,
            "n_vmr_cpgs_assayed": len(merged),
            "n_with_external_support": n_support,
            "fraction_with_external_support": n_support / len(merged) if len(merged) else 0.0,
            "assay_universe_complete": bool(merged["assay_universe_complete"].all()) if len(merged) else False,
            "n_vmrs_with_any_support": int(merged.loc[merged["external_meqtl_support"] == 1, "vmr_id"].nunique()),
            "harmonized_path": str(HARM / f"{resource_id}.hg38.tsv.gz"),
            "region_support_path": str(out),
        })
        print(f"  {region}: {n_support:,}/{len(merged):,} CpGs with external support")
    return pd.DataFrame(summaries)


def write_harmonized(resource_id: str, df: pd.DataFrame, summary: pd.DataFrame) -> None:
    HARM.mkdir(parents=True, exist_ok=True)
    out = HARM / f"{resource_id}.hg38.tsv.gz"
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    summary_path = HARM / f"{resource_id}.overlap_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote {out} ({len(df):,} rows)")
    print(f"Wrote {summary_path}")


def harmonize_jaffe(fdr: float, region_maps: dict[str, pd.DataFrame]) -> None:
    """Harmonize the Jaffe DLPFC 450K meQTL catalog against the 450K manifest.

    NOTE: `GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv.gz` is a *significant-results*
    table -- "allPairs" means every SNP partner of each significant CpG, not every pair
    tested (its maximum FDR is ~1e-27). It therefore cannot supply the tested universe;
    the 450K manifest does. See PHASE3_DIAGNOSIS.md.
    """
    path = RAW / "jaffe_dlpfc_450k_meqtl" / "geo" / "GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv.gz"
    if not path.exists():
        raise SystemExit(f"Missing Jaffe GEO file: {path}")
    print(f"Reading {path}")
    usecols = [
        "snps", "cpg", "statistic", "pvalue", "FDR", "beta",
        "snpChr", "snpPos", "snpRsNum", "methChr", "methPos",
    ]
    # Reduce each chunk to the best pair per CpG and merge incrementally to bound memory.
    lead = None
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=500_000):
        best = chunk.sort_values(["FDR", "pvalue"]).drop_duplicates("cpg", keep="first")
        lead = best if lead is None else pd.concat([lead, best], ignore_index=True)
        lead = lead.sort_values(["FDR", "pvalue"]).drop_duplicates("cpg", keep="first")
    if lead is None or lead.empty:
        raise SystemExit("Jaffe all-pairs table is empty")
    lead = lead.loc[lead["FDR"] <= fdr].copy()
    print(f"Jaffe significant CpGs at FDR<={fdr}: {lead['cpg'].nunique():,}")

    sig = pd.DataFrame({
        "probe_id": lead["cpg"].astype(str),
        "external_pvalue": lead["pvalue"],
        "external_fdr": lead["FDR"],
        "external_beta": lead["beta"],
        "lead_snp": lead["snpRsNum"].fillna(lead["snps"]),
        "snp_chrom": lead["snpChr"],
        "snp_pos": lead["snpPos"],
    })
    out, stats = build_array_resource("jaffe_dlpfc_450k_meqtl", "DLPFC", sig, fdr)
    summary = annotate_regions(out, "jaffe_dlpfc_450k_meqtl", region_maps)
    for k, v in stats.items():
        summary[k] = v
    write_harmonized("jaffe_dlpfc_450k_meqtl", out, summary)


def harmonize_schulz(fdr: float, region_maps: dict[str, pd.DataFrame]) -> None:
    path = RAW / "schulz_hippocampus_array_meqtl" / "nature_supplements" / "41467_2017_1818_MOESM4_ESM.xlsb"
    if not path.exists():
        raise SystemExit(f"Missing Schulz meQTL table: {path}")
    print(f"Reading {path}")
    df = pd.read_excel(path, sheet_name="meQTLs", engine="pyxlsb")
    # This supplement is a discovery catalog. It does not itself carry the tested-probe
    # universe, but Schulz is a 450K study (>99.9% of its probes are on the manifest),
    # so the manifest supplies the negatives. See PHASE3_DIAGNOSIS.md.
    df = df.loc[df["meQTL_FDR"] <= fdr].copy()
    print(f"Schulz significant CpGs at FDR<={fdr}: {df['CpG_meQTL_Probe'].nunique():,}")
    df = df.sort_values(["meQTL_FDR", "meQTL_P"], ascending=True)
    lead = df.drop_duplicates(subset=["CpG_meQTL_Probe"], keep="first").copy()

    sig = pd.DataFrame({
        "probe_id": lead["CpG_meQTL_Probe"].astype(str),
        "external_pvalue": lead["meQTL_P"],
        "external_fdr": lead["meQTL_FDR"],
        "external_beta": pd.NA,
        "lead_snp": lead["SNP_meQTL"],
        "snp_chrom": pd.NA,
        "snp_pos": lead["SNP_bp_hg19"],
    })
    out, stats = build_array_resource(
        "schulz_hippocampus_array_meqtl", "hippocampus", sig, fdr
    )
    summary = annotate_regions(out, "schulz_hippocampus_array_meqtl", region_maps)
    for k, v in stats.items():
        summary[k] = v
    write_harmonized("schulz_hippocampus_array_meqtl", out, summary)


def _read_brainseq_scz_table(path: Path, tissue: str) -> pd.DataFrame:
    raw = pd.read_excel(path, header=None)
    # Row 0 = title, row 1 = header names
    header = [str(x).strip() for x in raw.iloc[1].tolist()]
    df = raw.iloc[2:].copy()
    df.columns = header
    df = df.dropna(how="all")
    df["tissue_region"] = tissue
    return df


def harmonize_brainseq_scz_subset(fdr: float, region_maps: dict[str, pd.DataFrame]) -> None:
    """Interim BrainSeq resource from Nature SCZ-risk SNP meQTL supplements.

    Not a genome-wide cis catalog. Full results remain blocked on Synapse.
    """
    root = RAW / "brainseq_wgbs_meqtl" / "nature_supplements"
    tables = [
        (root / "41467_2021_25517_MOESM5_ESM.xlsx", "DLPFC"),
        (root / "41467_2021_25517_MOESM6_ESM.xlsx", "hippocampus"),
    ]
    frames = []
    for path, tissue in tables:
        if not path.exists():
            raise SystemExit(f"Missing BrainSeq supplement: {path}")
        frames.append(_read_brainseq_scz_table(path, tissue))
    df = pd.concat(frames, ignore_index=True)
    # cpg like chr14.103850655
    parts = df["cpg"].astype(str).str.extract(r"^(chr[\w]+)\.(\d+)$")
    df["chrom"] = parts[0]
    df["pos_1based"] = pd.to_numeric(parts[1], errors="coerce")
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["beta"] = pd.to_numeric(df["beta"], errors="coerce")
    df = df.dropna(subset=["chrom", "pos_1based"])
    df = df.loc[df["FDR"] <= fdr].copy()
    print(f"BrainSeq SCZ-subset pairs FDR<={fdr}: {len(df):,}; unique CpGs {df.groupby(['chrom','pos_1based']).ngroups:,}")

    df = df.sort_values(["FDR", "pvalue"], ascending=True)
    lead = df.drop_duplicates(subset=["chrom", "pos_1based", "tissue_region"], keep="first").copy()
    # Genome build: Perzel Mandell 2021 WGBS used GRCh38
    out = pd.DataFrame({
        "resource_id": "brainseq_wgbs_meqtl_scz_subset",
        "tissue_region": lead["tissue_region"],
        "phenotype_id": lead["chrom"].astype(str) + "_" + lead["pos_1based"].astype(int).astype(str),
        "chrom": lead["chrom"],
        "pos_1based": lead["pos_1based"].astype(int),
        "probe_id": lead["cpg"],
        "genome_build_source": "hg38",
        "genome_build_harmonized": "hg38",
        "lift_status": "native_hg38",
        "external_pvalue": lead["pvalue"],
        "external_fdr": lead["FDR"],
        "external_beta": lead["beta"],
        "lead_snp": lead["snps"],
        "snp_chrom": lead.get("snpChr"),
        "snp_pos": lead.get("snpPos"),
        "external_meqtl_support": 1,
        "external_assayed": 1,
        "assay_universe_complete": False,
        "notes": "SCZ-risk SNP meQTL supplement only; full cis catalog requires Synapse syn25992404",
    })
    summary = annotate_regions(out, "brainseq_wgbs_meqtl_scz_subset", region_maps)
    summary["catalog_scope"] = "scz_risk_snp_meqtl_supplement_only"
    summary["full_catalog_blocker"] = "SynapseAuthenticationError syn25992404"
    write_harmonized("brainseq_wgbs_meqtl_scz_subset", out, summary)

    # Also write a placeholder note for the primary resource id
    note = HARM / "brainseq_wgbs_meqtl.PENDING_SYNAPSE.txt"
    note.write_text(
        "Full BrainSeq WGBS cis-meQTL catalogs are on Synapse syn25992404 "
        "(PsychENCODE / DOI 10.7303/syn25992404).\n"
        "Current ~/.synapseConfig token is invalid (SynapseAuthenticationError).\n"
        "Interim harmonized file: brainseq_wgbs_meqtl_scz_subset.hg38.tsv.gz "
        "(Nature SCZ-risk SNP–CpG tables only; not genome-wide).\n"
        "After refreshing the Synapse token, re-download full results into "
        "raw/brainseq_wgbs_meqtl/synapse/ and extend this script.\n"
    )
    print(f"Wrote {note}")


def update_checklist(resources_done: list[str]) -> None:
    path = BASE / "download_checklist.tsv"
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    # Map subset / complete status
    status = {
        "jaffe_dlpfc_450k_meqtl": (
            True,
            str(HARM / "jaffe_dlpfc_450k_meqtl.hg38.tsv.gz"),
            "GEO GSE74193 significant meQTLs + hg38 harmonized; "
            "tested-negative universe from 450K manifest",
        ),
        "schulz_hippocampus_array_meqtl": (
            True,
            str(HARM / "schulz_hippocampus_array_meqtl.hg38.tsv.gz"),
            "Nature MOESM4 meQTLs + hg38 harmonized; "
            "tested-negative universe from 450K manifest",
        ),
        "brainseq_wgbs_meqtl": (
            False,
            str(HARM / "brainseq_wgbs_meqtl_scz_subset.hg38.tsv.gz"),
            "Nature SCZ-risk subset harmonized; full cis catalog blocked on Synapse token syn25992404",
        ),
    }
    rows = []
    for _, row in df.iterrows():
        rid = row["resource_id"]
        if rid in status:
            done, harm, notes = status[rid]
            row = row.copy()
            row["download_complete"] = str(done).lower()
            row["harmonized_path"] = harm
            row["notes"] = notes
        rows.append(row.to_dict())
    write_tsv(path, rows, list(df.columns))
    print(f"Updated {path}")


def main() -> None:
    args = parse_args()
    HARM.mkdir(parents=True, exist_ok=True)
    if not LIFTOVER.exists():
        raise SystemExit(f"liftOver not found: {LIFTOVER}")
    region_maps = load_cpg_vmr_maps(args.regions)
    if not region_maps:
        raise SystemExit("No region cpg_vmr_map tables found")

    handlers = {
        "jaffe_dlpfc_450k_meqtl": lambda: harmonize_jaffe(args.fdr, region_maps),
        "schulz_hippocampus_array_meqtl": lambda: harmonize_schulz(args.fdr, region_maps),
        "brainseq_wgbs_meqtl_scz_subset": lambda: harmonize_brainseq_scz_subset(args.fdr, region_maps),
        "brainseq_wgbs_meqtl": lambda: harmonize_brainseq_scz_subset(args.fdr, region_maps),
    }
    for rid in args.resources:
        if rid not in handlers:
            raise SystemExit(f"Unknown resource: {rid}")
        print(f"\n=== Harmonizing {rid} ===")
        handlers[rid]()
    update_checklist(args.resources)


if __name__ == "__main__":
    main()
