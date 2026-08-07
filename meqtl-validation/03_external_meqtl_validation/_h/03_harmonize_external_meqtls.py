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
    """Join external significant CpGs to each region's VMR CpG map."""
    sig = cpg_level.loc[cpg_level["external_meqtl_support"] == 1].copy()
    sig = sig.dropna(subset=["chrom", "pos_1based"])
    sig["chrom"] = sig["chrom"].astype(str)
    sig["pos_1based"] = sig["pos_1based"].astype(int)
    # Keep one row per CpG (strongest FDR)
    if "external_fdr" in sig.columns:
        sig = sig.sort_values("external_fdr", ascending=True)
    sig = sig.drop_duplicates(subset=["chrom", "pos_1based"], keep="first")

    summaries = []
    for region, mmap in region_maps.items():
        merged = mmap.merge(
            sig[["chrom", "pos_1based", "external_meqtl_support", "external_fdr", "external_pvalue", "lead_snp", "probe_id"]],
            on=["chrom", "pos_1based"],
            how="left",
        )
        merged["external_meqtl_support"] = merged["external_meqtl_support"].fillna(0).astype(int)
        merged["resource_id"] = resource_id
        merged["region"] = region
        out = HARM / f"{resource_id}.{region}.vmr_support.tsv.gz"
        merged.to_csv(out, sep="\t", index=False, compression="gzip")
        n_support = int(merged["external_meqtl_support"].sum())
        summaries.append({
            "resource_id": resource_id,
            "region": region,
            "n_vmr_cpgs": len(merged),
            "n_with_external_support": n_support,
            "fraction_with_external_support": n_support / len(merged) if len(merged) else 0.0,
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
    path = RAW / "jaffe_dlpfc_450k_meqtl" / "geo" / "GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv.gz"
    if not path.exists():
        raise SystemExit(f"Missing Jaffe GEO file: {path}")
    print(f"Reading {path}")
    usecols = [
        "snps", "cpg", "statistic", "pvalue", "FDR", "beta",
        "snpChr", "snpPos", "snpRsNum", "methChr", "methPos",
    ]
    # Stream filter to significant pairs to keep memory reasonable
    chunks = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=500_000):
        chunk = chunk.loc[chunk["FDR"] <= fdr]
        if not chunk.empty:
            chunks.append(chunk)
    if not chunks:
        raise SystemExit("No Jaffe pairs passed FDR filter")
    df = pd.concat(chunks, ignore_index=True)
    print(f"Jaffe pairs FDR<={fdr}: {len(df):,}; unique CpGs {df['cpg'].nunique():,}")

    # Lead SNP per CpG = lowest FDR (then lowest p)
    df = df.sort_values(["FDR", "pvalue"], ascending=True)
    lead = df.drop_duplicates(subset=["cpg"], keep="first").copy()
    lift = lift_positions(lead, "methChr", "methPos")
    lead = lead.merge(lift, on=["methChr", "methPos"], how="left")

    out = pd.DataFrame({
        "resource_id": "jaffe_dlpfc_450k_meqtl",
        "tissue_region": "DLPFC",
        "phenotype_id": lead.apply(
            lambda r: f"{r['chrom_hg38']}_{int(r['pos_hg38'])}" if pd.notna(r["pos_hg38"]) else "",
            axis=1,
        ),
        "chrom": lead["chrom_hg38"],
        "pos_1based": lead["pos_hg38"],
        "probe_id": lead["cpg"],
        "genome_build_source": "hg19",
        "genome_build_harmonized": "hg38",
        "lift_status": lead["lift_status"],
        "external_pvalue": lead["pvalue"],
        "external_fdr": lead["FDR"],
        "external_beta": lead["beta"],
        "lead_snp": lead["snpRsNum"].fillna(lead["snps"]),
        "snp_chrom": lead["snpChr"],
        "snp_pos": lead["snpPos"],
        "external_meqtl_support": (lead["lift_status"] == "mapped").astype(int),
    })
    # Only mapped CpGs count as support
    out.loc[out["lift_status"] != "mapped", "external_meqtl_support"] = 0
    summary = annotate_regions(out, "jaffe_dlpfc_450k_meqtl", region_maps)
    # Add lift stats to summary
    lift_stats = {
        "n_unique_cpgs_source": int(lead["cpg"].nunique()),
        "n_lift_mapped": int((lead["lift_status"] == "mapped").sum()),
        "n_lift_unmapped": int((lead["lift_status"] != "mapped").sum()),
        "fdr_threshold": fdr,
    }
    for k, v in lift_stats.items():
        summary[k] = v
    write_harmonized("jaffe_dlpfc_450k_meqtl", out, summary)


def harmonize_schulz(fdr: float, region_maps: dict[str, pd.DataFrame]) -> None:
    path = RAW / "schulz_hippocampus_array_meqtl" / "nature_supplements" / "41467_2017_1818_MOESM4_ESM.xlsb"
    if not path.exists():
        raise SystemExit(f"Missing Schulz meQTL table: {path}")
    print(f"Reading {path}")
    df = pd.read_excel(path, sheet_name="meQTLs", engine="pyxlsb")
    df = df.loc[df["meQTL_FDR"] <= fdr].copy()
    print(f"Schulz pairs FDR<={fdr}: {len(df):,}; unique CpGs {df['CpG_meQTL_Probe'].nunique():,}")
    df = df.sort_values(["meQTL_FDR", "meQTL_P"], ascending=True)
    lead = df.drop_duplicates(subset=["CpG_meQTL_Probe"], keep="first").copy()
    lift = lift_positions(lead, "CpG_chr_hg19", "CpG_bp_hg19")
    lead = lead.merge(lift, on=["CpG_chr_hg19", "CpG_bp_hg19"], how="left")

    out = pd.DataFrame({
        "resource_id": "schulz_hippocampus_array_meqtl",
        "tissue_region": "hippocampus",
        "phenotype_id": lead.apply(
            lambda r: f"{r['chrom_hg38']}_{int(r['pos_hg38'])}" if pd.notna(r["pos_hg38"]) else "",
            axis=1,
        ),
        "chrom": lead["chrom_hg38"],
        "pos_1based": lead["pos_hg38"],
        "probe_id": lead["CpG_meQTL_Probe"],
        "genome_build_source": "hg19",
        "genome_build_harmonized": "hg38",
        "lift_status": lead["lift_status"],
        "external_pvalue": lead["meQTL_P"],
        "external_fdr": lead["meQTL_FDR"],
        "external_beta": pd.NA,
        "lead_snp": lead["SNP_meQTL"],
        "snp_chrom": pd.NA,
        "snp_pos": lead["SNP_bp_hg19"],
        "external_meqtl_support": (lead["lift_status"] == "mapped").astype(int),
    })
    out.loc[out["lift_status"] != "mapped", "external_meqtl_support"] = 0
    summary = annotate_regions(out, "schulz_hippocampus_array_meqtl", region_maps)
    for k, v in {
        "n_unique_cpgs_source": int(lead["CpG_meQTL_Probe"].nunique()),
        "n_lift_mapped": int((lead["lift_status"] == "mapped").sum()),
        "n_lift_unmapped": int((lead["lift_status"] != "mapped").sum()),
        "fdr_threshold": fdr,
    }.items():
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
            "GEO GSE74193 allPairs downloaded + hg38 harmonized",
        ),
        "schulz_hippocampus_array_meqtl": (
            True,
            str(HARM / "schulz_hippocampus_array_meqtl.hg38.tsv.gz"),
            "Nature MOESM4 meQTLs downloaded + hg38 harmonized",
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
