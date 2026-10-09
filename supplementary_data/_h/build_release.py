#!/usr/bin/env python3
"""Tier every file of every accepted run for release: git, LFS, Zenodo, or none.

The rule, set when the results were first shared with collaborators:

  * Every file the manuscript is written from is pushed to GitHub -- plain git
    when it is 15 MB or smaller, Git LFS above that. Zenodo carries a mirror
    copy for the DOI and is never the only home of such a file.
  * Reproducibility extras that no manuscript text draws on (run provenance,
    per-VMR tables nobody cites, annotations, LD scores) go to Zenodo only,
    one zip per module.
  * Individual-level data -- methylation matrices, genotypes, covariates,
    phenotype tables -- goes nowhere. It is governed by dbGaP phs000979.
    Bulk intermediates (per-task shards, nominal meQTL output, coloc region
    shards, copies of public GWAS) stay on Quest; they are regenerable from
    `_h/` and no claim rests on them.

A "manuscript file" is any of:
  1. a `source_table` resolved from the newest accepted Module 11 run's
     `tables/manuscript-number-registry.tsv`;
  2. a deposit-ready row of `supplementary_data/supplementary_data_manifest.tsv`;
  3. an accepted run's `manifest.tsv` and its decision, gate, claim and
     interpretation tables, so every number can be tied to a run and commit;
  4. the whole accepted Module 11 run (figures, source data, tables).
The individual-level rule is checked first and wins over all four.

Accepted runs come from `00_shared/gates.R` (`read_accepted_runs()` and
`read_accepted_nongating_runs()`), the parser of record, through Rscript.

  python3 supplementary_data/_h/build_release.py            # plan: write the manifests
  python3 supplementary_data/_h/build_release.py --zip      # also build zenodo/{module}.zip
  python3 supplementary_data/_h/build_release.py --stage    # git add -f the git/LFS tiers
  python3 supplementary_data/_h/build_release.py --check    # verify, exit non-zero on a problem
"""
import argparse, csv, glob, gzip, hashlib, io, os, re, subprocess, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from verify_manifest import expand_braces  # noqa: E402

GIT_MAX = 15 * 1024 * 1024          # above this a manuscript file goes to LFS
ZENODO_FILE_MAX = 1024 ** 3         # a single file above this stays on Quest
RELEASE_MANIFEST = os.path.join(ROOT, "supplementary_data", "release_manifest.tsv")
ZENODO_DIR = os.path.join(ROOT, "zenodo")
ZENODO_MANIFEST = os.path.join(ZENODO_DIR, "zenodo_manifest.tsv")
GITATTR = os.path.join(ROOT, ".gitattributes")
LFS_BEGIN = "# >>> release LFS tier (supplementary_data/_h/build_release.py) >>>"
LFS_END = "# <<< release LFS tier <<<"
R_BIN = os.environ.get("V2_RSCRIPT", "/projects/p32505/opt/envs/epigenomics/bin/Rscript")
REGIONS = ("caudate", "dlpfc", "hippocampus")

# ---------------------------------------------------------------- the rules
# Individual-level data, by path segment or extension. Checked first.
INDIVIDUAL_SEGMENTS = {"cpg", "plink_format", "covs", "phenotypes", "geno_stage",
                       "tested_meth", "pca"}
# Donor lists that pair a BrNum with its array barcode. The barcode is stripped
# before any release (PI decision 2026-09-20), so these files stay on Quest.
INDIVIDUAL_NAMES = {"donors_plink.txt"}
INDIVIDUAL_EXT = (".pgen", ".pvar", ".psam", ".bim", ".fam", ".phen", ".covar",
                  ".qcovar", ".grm.bin", ".grm.N.bin", ".grm.id", ".bgen", ".vcf",
                  ".vcf.gz")
# Bulk intermediates and external copies: Quest only, whatever their tier.
BULK_SEGMENTS = {"task_rows", "oof", "nominal", "subsample", "regions", "sumstats",
                 "gwas", "checkpoint", "work", "status", "logs", "code", "inputs",
                 "ld", "raw"}
# What the Zenodo-only tier may contain, by first path segment under the run.
ZENODO_SEGMENTS = {"results", "summary", "annotation", "ldscores", "links",
                   "excluded", "figures", "vmr", "qc", "config", "coloc",
                   "combined", "provenance"}
ZENODO_EXT = (".tsv", ".tsv.gz", ".txt", ".json", ".results", ".log", ".bed",
              ".gz", ".png", ".pdf", ".M", ".M_5_50")
DECISION_RE = re.compile(r"(decision|gate|claims|interpretation|reading|"
                         r"acceptance|reconciliation|qc|overall|by-h2|stratified)",
                         re.I)
DONOR_COL_RE = re.compile(r"^Br\d+$")


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)


def accepted_runs():
    """[(module, run_id, gating)] from the R parser of record."""
    code = r'''
suppressMessages(source("00_shared/load.R"))
mods <- sort(list.dirs(".", recursive = FALSE, full.names = FALSE))
mods <- mods[grepl("^[0-9][0-9]", mods) & file.exists(file.path(mods, "README.md"))]
for (m in mods) {
  for (g in c(TRUE, FALSE)) {
    d <- tryCatch(if (g) read_accepted_runs(m) else read_accepted_nongating_runs(m),
                  error = function(e) NULL)
    if (!is.null(d) && nrow(d)) for (r in d$run_id) cat(m, r, g, sep = "\t", "\n")
  }
}'''
    r = subprocess.run([R_BIN, "-e", code], capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        sys.exit("read_accepted_runs() failed:\n" + r.stderr[-2000:])
    out = []
    for line in r.stdout.splitlines():
        p = line.strip().split("\t")
        if len(p) == 3 and os.path.isdir(os.path.join(ROOT, p[0], "_m", "runs", p[1])):
            out.append((p[0], p[1], p[2] == "TRUE"))
    return out


def decision_runs(accepted_ids):
    """Runs a locked config pins as decision evidence (non-comment YAML values)."""
    tok = re.compile(r"[a-z0-9]+(?:-[A-Za-z0-9_.]+)+-\d{8}(?:-?[a-z])?(?![A-Za-z0-9])")
    found = set()
    for f in glob.glob(os.path.join(ROOT, "config", "*.yml")):
        for ln in open(f, encoding="utf-8", errors="replace"):
            if ln.lstrip().startswith("#"):
                continue
            found |= set(tok.findall(re.sub(r"\s#.*$", "", ln)))
    out = []
    for rid in sorted(found - accepted_ids):
        for d in glob.glob(os.path.join(ROOT, "*", "_m", "runs", rid)):
            out.append((os.path.relpath(d, ROOT).split(os.sep)[0], rid, False))
    return out


def walk(run_dir, pruned):
    for dp, dn, fn in os.walk(run_dir):
        rel = os.path.relpath(dp, run_dir)
        segs = [] if rel == "." else rel.split(os.sep)
        # Never descend into a bulk or individual-level tree: they hold most of
        # the 1.5 M files, and nothing inside them can be released. Each pruned
        # directory is still recorded, once, so the exclusion is auditable.
        for d in dn:
            if d in INDIVIDUAL_SEGMENTS:
                pruned.append((os.path.join(dp, d), "individual-level data (dbGaP phs000979)"))
            elif d in BULK_SEGMENTS:
                pruned.append((os.path.join(dp, d), "bulk or regenerable intermediate (Quest)"))
        dn[:] = [d for d in dn if d not in INDIVIDUAL_SEGMENTS | BULK_SEGMENTS]
        for f in fn:
            yield os.path.join(dp, f), segs


def is_individual(path, segs):
    if INDIVIDUAL_SEGMENTS & set(segs) or path.endswith(INDIVIDUAL_EXT) or \
            os.path.basename(path) in INDIVIDUAL_NAMES:
        return True
    if path.endswith(".bed") and os.path.exists(path[:-4] + ".bim"):
        return True  # a PLINK genotype .bed, not an interval BED
    if path.endswith((".tsv", ".tsv.gz", ".txt", ".csv")):
        try:
            op = gzip.open if path.endswith(".gz") else open
            with op(path, "rt", errors="replace") as fh:
                head = fh.readline().rstrip("\n").split("\t")
            # A donor-by-feature matrix: many columns that are donor IDs.
            if sum(bool(DONOR_COL_RE.match(c)) for c in head) >= 20:
                return True
        except (OSError, EOFError, UnicodeDecodeError):
            pass
    return False


def registry_files(fig_run):
    """Files resolved from the Module 11 number registry's source_table strings."""
    reg = os.path.join(ROOT, "11_integrated_manuscript_outputs", "_m", "runs", fig_run,
                       "tables", "manuscript-number-registry.tsv")
    out = set()
    if not os.path.exists(reg):
        return out
    with open(reg) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        runs_ = [x for x in (r.get("source_run_id") or "").split(";") if x]
        for tok in (r.get("source_table") or "").split(" + "):
            tok = re.sub(r"^\d\d[a-z]?: ", "", tok.strip())
            tok = re.sub(r"^\d\d[a-z]?_[a-z_]+ ", "", tok)
            tok = tok.split(" ")[0]
            for rid in runs_:
                for rd in glob.glob(os.path.join(ROOT, "*", "_m", "runs", rid)):
                    mod_root = os.path.dirname(os.path.dirname(os.path.dirname(rd)))
                    for pat in expand_braces(tok.replace("{region}", "{" + ",".join(REGIONS) + "}")
                                             .replace("{cohort}", "AA")):
                        for base in (rd, os.path.join(rd, "results"), mod_root):
                            p = os.path.normpath(os.path.join(base, pat))
                            if os.path.isfile(p):
                                out.add(os.path.relpath(p, ROOT))
    return out


def sd_files():
    """Deposit-ready Supplementary Data locations, with their SD numbers."""
    man = os.path.join(ROOT, "supplementary_data", "supplementary_data_manifest.tsv")
    ok = {"ready", "needs_deid", "pending_acceptance", "transform_on_deposit"}
    out = {}
    with open(man) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] not in ok:
                continue
            for pat in expand_braces(r["location"]):
                for p in glob.glob(os.path.join(ROOT, pat)):
                    if os.path.isfile(p):
                        out.setdefault(os.path.relpath(p, ROOT), set()).add(r["sd"])
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def plan():
    acc = accepted_runs()
    acc_ids = {r for _, r, _ in acc}
    runs_ = [(m, r, g, "accepted") for m, r, g in acc] + \
            [(m, r, g, "decision_evidence") for m, r, g in decision_runs(acc_ids)]
    figs = sorted(r for m, r, _ in acc if m == "11_integrated_manuscript_outputs")
    fig_run = figs[-1] if figs else None
    reg = registry_files(fig_run) if fig_run else set()
    sd = sd_files()

    rows = []
    # The tracked _m/combined/ deliverables are manuscript files already in git;
    # they are listed so the Zenodo zips mirror them too.
    for rel in run(["git", "ls-files", "--", "*/_m/combined/*"]).stdout.split("\n"):
        if not rel or not os.path.isfile(os.path.join(ROOT, rel)):
            continue
        size = os.path.getsize(os.path.join(ROOT, rel))
        rows.append({"module": rel.split("/")[0], "run_id": "combined",
                     "run_kind": "combined", "path": rel, "bytes": str(size),
                     "tier": "git" if size <= GIT_MAX else "lfs",
                     "reason": "cross-region deliverable (_m/combined)",
                     "sd": ",".join(sorted(sd.get(rel, ()), key=int))})
    # Supplementary Data that lives outside module runs (inputs/cell_proportions)
    # is tracked already; list it so the Zenodo mirror carries it.
    listed = {r["path"] for r in rows}
    tracked_all = set(run(["git", "ls-files"]).stdout.split("\n"))
    for rel in sorted(sd):
        if "/_m/runs/" in rel or rel in listed or rel not in tracked_all:
            continue
        size = os.path.getsize(os.path.join(ROOT, rel))
        rows.append({"module": rel.split("/")[0], "run_id": "tracked",
                     "run_kind": "tracked", "path": rel, "bytes": str(size),
                     "tier": "git" if size <= GIT_MAX else "lfs",
                     "reason": "Supplementary Data " + ",".join(sorted(sd[rel], key=int)),
                     "sd": ",".join(sorted(sd[rel], key=int))})
    for mod, rid, gating, kind in runs_:
        rd = os.path.join(ROOT, mod, "_m", "runs", rid)
        pruned = []
        for p, segs in walk(rd, pruned):
            rel = os.path.relpath(p, ROOT)
            name = os.path.basename(p)
            size = os.path.getsize(p)
            top = segs[0] if segs else ""
            manuscript, why = False, ""
            if rel in reg:
                manuscript, why = True, "figure source table (Module 11 registry)"
            elif rel in sd:
                manuscript, why = True, "Supplementary Data " + ",".join(sorted(sd[rel], key=int))
            elif mod == "11_integrated_manuscript_outputs" and rid == fig_run and \
                    (top in {"figures", "source_data", "tables"} or not segs):
                manuscript, why = True, "accepted figure run"
            elif not segs and name == "manifest.tsv":
                manuscript, why = True, "run manifest"
            elif top in {"", "results", "combined"} and len(segs) <= 1 \
                    and DECISION_RE.search(name) and size <= GIT_MAX:
                manuscript, why = True, "decision/gate table"

            if is_individual(p, segs):
                tier, why = "none", "individual-level data (dbGaP phs000979)"
            elif name.endswith(".svg"):
                tier, why = "none", "SVG duplicate of a PDF/PNG figure"
            elif manuscript:
                tier = "git" if size <= GIT_MAX else "lfs"
            elif (top in ZENODO_SEGMENTS or not segs) and name.endswith(ZENODO_EXT) \
                    and size <= ZENODO_FILE_MAX:
                tier, why = "zenodo", "reproducibility extra (no manuscript text uses it)"
            else:
                tier, why = "none", "bulk or regenerable intermediate (Quest)"
            rows.append({"module": mod, "run_id": rid, "run_kind": kind,
                         "path": rel, "bytes": str(size), "tier": tier, "reason": why,
                         "sd": ",".join(sorted(sd.get(rel, ()), key=int))})
        for d, why in pruned:
            rows.append({"module": mod, "run_id": rid, "run_kind": kind,
                         "path": os.path.relpath(d, ROOT) + "/", "bytes": "0",
                         "tier": "none", "reason": why + "; directory not listed",
                         "sd": ""})
    return rows, fig_run


def write_tsv(path, rows, cols):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def write_lfs_rules(rows):
    lfs = sorted(r["path"] for r in rows if r["tier"] == "lfs")
    text = open(GITATTR).read() if os.path.exists(GITATTR) else ""
    text = re.sub(re.escape(LFS_BEGIN) + r".*?" + re.escape(LFS_END) + r"\n?", "", text,
                  flags=re.S)
    block = "\n".join([LFS_BEGIN] + [f"{p.replace(' ', '[[:space:]]')} filter=lfs "
                                     f"diff=lfs merge=lfs -text" for p in lfs] + [LFS_END])
    with open(GITATTR, "w") as fh:
        fh.write(text.rstrip("\n") + "\n" + block + "\n")


def build_zips(rows):
    os.makedirs(ZENODO_DIR, exist_ok=True)
    by_mod = {}
    for r in rows:
        if r["tier"] in ("git", "lfs", "zenodo"):  # Zenodo mirrors the manuscript files
            by_mod.setdefault(r["module"], []).append(r)
    out = []
    for mod, rs in sorted(by_mod.items()):
        z = os.path.join(ZENODO_DIR, f"{mod}.zip")
        tmp = z + ".tmp"
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as zf:
            for r in sorted(rs, key=lambda x: x["path"]):
                zf.write(os.path.join(ROOT, r["path"]), r["path"])
        os.replace(tmp, z)
        sds = sorted({s for r in rs for s in r["sd"].split(",") if s}, key=int)
        out.append({"archive": os.path.basename(z), "module": mod,
                    "bytes": str(os.path.getsize(z)), "sha256": sha256(z),
                    "n_members": str(len(rs)),
                    "n_manuscript_members": str(sum(r["tier"] != "zenodo" for r in rs)),
                    "sd_numbers": ",".join(sds),
                    "runs": ";".join(sorted({r["run_id"] for r in rs}))})
        print(f"[zip] {os.path.basename(z)}: {len(rs)} files, "
              f"{os.path.getsize(z) / 1e6:.1f} MB")
    write_tsv(ZENODO_MANIFEST, out, list(out[0].keys()) if out else ["archive"])


def stage(rows):
    paths = [r["path"] for r in rows if r["tier"] in ("git", "lfs")]
    for i in range(0, len(paths), 500):
        r = run(["git", "add", "-f", "--"] + paths[i:i + 500])
        if r.returncode != 0:
            sys.exit(r.stderr)
    run(["git", "add", ".gitattributes", os.path.relpath(RELEASE_MANIFEST, ROOT)])
    print(f"[stage] {len(paths)} files staged")


def check(rows):
    bad = []
    tracked = set(run(["git", "ls-files"]).stdout.split("\n"))
    lfs = set(l.split(" ")[-1] for l in run(["git", "lfs", "ls-files", "-n"]).stdout.split("\n") if l)
    for r in rows:
        p = r["path"]
        if r["tier"] in ("git", "lfs") and p not in tracked:
            bad.append(f"manuscript file not tracked: {p}")
        if r["tier"] == "lfs" and p not in lfs:
            bad.append(f"LFS-tier file not in LFS: {p}")
        if r["tier"] == "none" and p in tracked:
            bad.append(f"none-tier file is tracked: {p} ({r['reason']})")
    if os.path.exists(ZENODO_MANIFEST):
        with open(ZENODO_MANIFEST) as fh:
            for z in csv.DictReader(fh, delimiter="\t"):
                f = os.path.join(ZENODO_DIR, z["archive"])
                if not os.path.exists(f):
                    bad.append(f"zip missing: {z['archive']}")
                    continue
                members = set(zipfile.ZipFile(f).namelist())
                none_in = [r["path"] for r in rows if r["tier"] == "none" and r["path"] in members]
                if none_in:
                    bad.append(f"{z['archive']} contains {len(none_in)} none-tier files")
                if str(len(members)) != z["n_members"]:
                    bad.append(f"{z['archive']} member count differs from the manifest")
    for b in bad[:50]:
        print("PROBLEM", b)
    print(f"[check] {len(bad)} problem(s)")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", action="store_true")
    ap.add_argument("--stage", action="store_true")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    rows, fig_run = plan()
    cols = ["module", "run_id", "run_kind", "path", "bytes", "tier", "reason", "sd"]
    write_tsv(RELEASE_MANIFEST, rows, cols)
    tot = {}
    for r in rows:
        t = tot.setdefault(r["tier"], [0, 0])
        t[0] += 1
        t[1] += int(r["bytes"])
    print(f"figure run of record: {fig_run}")
    for t in ("git", "lfs", "zenodo", "none"):
        n, b = tot.get(t, [0, 0])
        print(f"  {t:<7} {n:8d} files {b / 1e9:9.3f} GB")
    write_lfs_rules(rows)
    if a.zip:
        build_zips(rows)
    if a.stage:
        stage(rows)
    if a.check:
        return check(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
