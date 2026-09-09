#!/usr/bin/env python3
"""Risk-variant -> CpG meQTL evidence for one chromosome (09).

Usage (inside an array task):
    python _h/03_risk_variant_cpg_meqtl.py --run-id scz-AA-caudate-YYYYMMDD --chrom 22

AGENTS.md 7.8: "link risk variants to corrected VMRs through explicit CpG meQTL
evidence". The evidence already exists -- the accepted Module 05 run wrote the
full nominal cis pass, every tested CpG x SNP pair, under
results/meqtl/nominal/. This stage does not remap anything. It selects the
prespecified risk variants, intersects them with the CpGs of the VMRs that
02_link_loci_to_vmrs.R linked to a locus, and writes that slice.

Risk variants are the published PGC3 index SNPs plus their LD proxies at
r2 >= config ld_r2_proxy_min, per config/analysis_thresholds.yml:phase7_scz.
LD is computed in the analysis cohort's own genotypes: the meQTL association
being annotated was estimated in those donors, so a proxy defined in a
different panel would not be a proxy for this test.

Multiple testing is NOT applied here. The FDR family is
`scz_risk_variant_cpg_pairs_per_region` -- it spans all autosomes of one region,
so correcting inside a per-chromosome array task would create 22 unrelated
families whose membership depended on how the work was parallelised.
03b_combine_risk_variant_tests.R pools the chromosomes and corrects once.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds
import pyarrow.compute as pc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_config import load_run_config  # noqa: E402

PLINK2 = "/projects/p32505/opt/bin/plink2"


def repo_root() -> Path:
    d = Path(__file__).resolve()
    while d != d.parent:
        if (d / ".git").is_dir():
            return d
        d = d.parent
    raise SystemExit("Could not locate repository root")


def read_manifest(run_dir: Path) -> dict:
    m = pd.read_csv(run_dir / "manifest.tsv", sep="\t", dtype=str)
    return dict(zip(m["field"], m["value"]))


def rsid_of(variant_id: str) -> str:
    """Module 05 variant IDs are chr_pos_ref_alt_rsID; the rsID is the tail."""
    tail = str(variant_id).rsplit("_", 1)[-1]
    return tail if tail.startswith("rs") else ""


def ld_proxies(index_snps: pd.DataFrame, chrom_label: str, run_dir: Path,
               root: Path, cohort: str, r2_min: float, window_bp: int,
               threads: int) -> pd.DataFrame:
    """Variants in LD with an index SNP at r2 >= r2_min, in the cohort panel.

    Returns a frame of index_snp / proxy rsID / r2. An index SNP absent from the
    panel yields no proxies and is recorded as such by the caller -- it is not
    an error, it is a coverage fact about the genotype data.
    """
    import yaml

    paths = yaml.safe_load((root / "config" / "paths.yml").read_text())
    arm = paths["genotype"]["AA" if cohort == "AA" else "all_individuals"]
    src_prefix = str(root / arm["pgen"])[: -len(".pgen")]

    work = run_dir / "coloc" / "ld" / chrom_label
    work.mkdir(parents=True, exist_ok=True)

    # plink2 refuses the distributed headerless .psam, so stage a headered copy
    # beside symlinked .pgen/.pvar rather than touching the shared source file.
    staged = work / "src"
    for ext in ("pgen", "pvar"):
        link = Path(f"{staged}.{ext}")
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(f"{src_prefix}.{ext}")
    psam_lines = ["#FID\tIID\tSEX"]
    for line in Path(f"{src_prefix}.psam").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        psam_lines.append(f"{f[0]}\t{f[1]}\t{f[2] if len(f) > 2 else 'NA'}")
    Path(f"{staged}.psam").write_text("\n".join(psam_lines) + "\n")

    snp_list = work / "index_snps.txt"
    snp_list.write_text("\n".join(index_snps["panel_variant_id"]) + "\n")

    out_prefix = work / "proxies"
    cmd = [
        PLINK2, "--pfile", str(staged), "--chr", chrom_label.replace("chr", ""),
        "--ld-snp-list", str(snp_list),
        "--r2-unphased", "--ld-window-r2", str(r2_min),
        "--ld-window-kb", str(int(window_bp / 1000)),
        "--threads", str(threads), "--out", str(out_prefix),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    vcor = Path(str(out_prefix) + ".vcor")
    if not vcor.exists():
        alt = list(work.glob("proxies.vcor*"))
        if not alt:
            raise SystemExit(
                f"plink2 produced no LD table for {chrom_label}:\n{res.stderr}")
        vcor = alt[0]
    ld = pd.read_csv(vcor, sep=r"\s+")
    cols = {c.upper(): c for c in ld.columns}
    id_a = cols.get("ID_A", cols.get("SNP_A"))
    id_b = cols.get("ID_B", cols.get("SNP_B"))
    r2c = cols.get("UNPHASED_R2", cols.get("R2", cols.get("PHASED_R2")))
    if not (id_a and id_b and r2c):
        raise SystemExit(f"Unexpected plink2 LD columns: {list(ld.columns)}")
    return pd.DataFrame({
        "panel_variant_id": ld[id_a].astype(str),
        "proxy_variant_id": ld[id_b].astype(str),
        "ld_r2": pd.to_numeric(ld[r2c], errors="coerce"),
    })


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    root = repo_root()
    run_dir = root / "09_schizophrenia_risk_application" / "_m" / "runs" / args.run_id
    if not run_dir.is_dir():
        raise SystemExit(f"No such run: {run_dir}")
    man = read_manifest(run_dir)
    cohort, region = man["cohort"], man["region"]
    cfg = load_run_config("schizophrenia", run_dir, root)

    chrom_label = f"chr{str(args.chrom).replace('chr', '')}"
    out_dir = run_dir / "results" / "risk_variant_tests"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_f = out_dir / f"{chrom_label}.tsv.gz"

    # ------------------------------------------------ the linked VMRs and CpGs
    links = pd.read_csv(run_dir / "results" / "locus-vmr-links.tsv", sep="\t")
    links = links[links["chrom"].astype(str) == chrom_label]
    index = pd.read_csv(run_dir / "results" / "scz-index-snps.tsv", sep="\t")
    index = index[(index["chrom"].astype(str) == chrom_label)
                  & index["locus_id"].notna()]

    if links.empty or index.empty:
        # A chromosome with no linked VMR or no assigned index SNP is a real,
        # recordable outcome. Write the empty slice so reconciliation sees the
        # task completed rather than inferring a silent failure.
        pd.DataFrame(columns=[
            "locus_id", "index_snp", "risk_variant_id", "ld_r2", "cpg_id",
            "vmr_id", "pval_nominal", "slope", "slope_se", "af",
            "start_distance", "chrom", "region", "population",
        ]).to_csv(out_f, sep="\t", index=False)
        print(f"{chrom_label}: no linked VMR or index SNP; wrote empty slice")
        return

    burden_run = man["upstream_cpg_meqtl_burden_run_id"]
    burden_dir = root / "05_cpg_meqtl_burden" / "_m" / "runs" / burden_run
    member = pd.read_csv(burden_dir / "results" / "tested-cpg-membership.tsv",
                         sep="\t")
    member = member[member["vmr_id"].isin(set(links["vmr_id"]))]
    if member.empty:
        pd.DataFrame(columns=["locus_id"]).to_csv(out_f, sep="\t", index=False)
        print(f"{chrom_label}: no tested CpG in any linked VMR")
        return
    cpg_ids = set(member["cpg_id"].astype(str))

    # ------------------------------------------------------- risk variant set
    nominal_f = (burden_dir / "results" / "meqtl" / "nominal"
                 / f"{chrom_label}.cis_qtl_pairs.{chrom_label}.parquet")
    if not nominal_f.is_file():
        raise SystemExit(
            f"Missing Module 05 nominal meQTL pairs for {chrom_label}: {nominal_f}. "
            "Module 09 does not remap meQTL; it consumes the accepted run.")

    # The panel variant IDs carry the rsID as their last underscore field, so
    # the published index rsIDs can be resolved to panel IDs without a
    # coordinate join (and therefore without a build assumption).
    pvar_ids = pd.read_parquet(nominal_f, columns=["variant_id"])
    pvar_ids = pd.DataFrame({"variant_id": pvar_ids["variant_id"].unique()})
    pvar_ids["rsid"] = pvar_ids["variant_id"].map(rsid_of)
    rs_to_panel = (pvar_ids[pvar_ids["rsid"].ne("")]
                   .drop_duplicates("rsid", keep="first")
                   .set_index("rsid")["variant_id"])

    index = index.copy()
    index["panel_variant_id"] = index["index_snp"].map(rs_to_panel)
    in_panel = index[index["panel_variant_id"].notna()]

    proxies = pd.DataFrame(columns=["panel_variant_id", "proxy_variant_id", "ld_r2"])
    if not in_panel.empty:
        proxies = ld_proxies(
            in_panel, chrom_label, run_dir, root, cohort,
            r2_min=float(cfg["testing"]["ld_r2_proxy_min"]),
            window_bp=int(cfg["testing"]["ld_proxy_window_bp"]),
            threads=args.threads)

    # The index SNP is its own proxy at r2 = 1, so a locus is never lost merely
    # because plink reported no partner above the threshold.
    self_rows = pd.DataFrame({
        "panel_variant_id": in_panel["panel_variant_id"],
        "proxy_variant_id": in_panel["panel_variant_id"],
        "ld_r2": 1.0,
    })
    proxies = (pd.concat([self_rows, proxies], ignore_index=True)
               .dropna(subset=["proxy_variant_id"])
               .sort_values("ld_r2", ascending=False)
               .drop_duplicates(["panel_variant_id", "proxy_variant_id"]))
    risk = in_panel[["locus_id", "index_snp", "panel_variant_id"]].merge(
        proxies, on="panel_variant_id", how="inner")
    risk = risk.rename(columns={"proxy_variant_id": "risk_variant_id"})
    if risk.empty:
        pd.DataFrame(columns=["locus_id"]).to_csv(out_f, sep="\t", index=False)
        print(f"{chrom_label}: no risk variant present in the genotype panel")
        return

    # --------------------------------------------------- the meQTL pair slice
    # Push both filters into the parquet scan. The nominal file is ~0.5-0.7 GB
    # per chromosome, and materialising it before filtering is what makes this
    # stage run out of memory.
    want_variants = set(risk["risk_variant_id"])
    scan = ds.dataset(nominal_f, format="parquet")
    tbl = scan.to_table(filter=(
        pc.field("variant_id").isin(list(want_variants))
        & pc.field("phenotype_id").isin(list(cpg_ids))))
    pairs = tbl.to_pandas()

    if pairs.empty:
        pd.DataFrame(columns=["locus_id"]).to_csv(out_f, sep="\t", index=False)
        print(f"{chrom_label}: risk variants present but no tested CpG pair")
        return

    pairs = pairs.rename(columns={"phenotype_id": "cpg_id",
                                  "variant_id": "risk_variant_id"})
    pairs = pairs.merge(risk, on="risk_variant_id", how="inner")
    pairs = pairs.merge(member[["cpg_id", "vmr_id"]].drop_duplicates(),
                        on="cpg_id", how="inner")
    # A CpG can sit in only one VMR, and a risk variant can proxy for only one
    # index SNP per locus; keep the strongest LD when a variant proxies several.
    pairs = (pairs.sort_values(["locus_id", "cpg_id", "risk_variant_id", "ld_r2"],
                               ascending=[True, True, True, False])
             .drop_duplicates(["locus_id", "cpg_id", "risk_variant_id"]))
    pairs["chrom"] = chrom_label
    pairs["region"] = region
    pairs["population"] = cohort

    keep = ["locus_id", "index_snp", "risk_variant_id", "ld_r2", "cpg_id",
            "vmr_id", "pval_nominal", "slope", "slope_se", "af",
            "start_distance", "chrom", "region", "population"]
    pairs[keep].to_csv(out_f, sep="\t", index=False)
    print(f"{chrom_label}: {len(pairs)} risk-variant x CpG pairs across "
          f"{pairs['locus_id'].nunique()} loci")


if __name__ == "__main__":
    main()
