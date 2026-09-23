#!/usr/bin/env python3
"""Mark superseded run directories for a later cleanup stage.

AGENTS.md 5.2 makes `_m/runs/{RUN_ID}/` immutable and 3 requires superseded
outputs to be removed from the submission branch once every downstream consumer
has a validated replacement. Between those two rules sits a directory tree that
nobody may edit and that nothing records as disposable. This script writes that
record: `DEPRECATED_RUNS.tsv`, one row per run directory that no accepted result
depends on.

It marks. It does not delete, and it must not learn to -- removal is a separate,
reviewed cleanup stage, because the cost of deleting a run that turns out to be
load-bearing is unrecoverable while the cost of keeping one is disk.

A run is PROTECTED when any of these holds:

  1. it appears in its module README's `## Accepted runs` table;
  2. it appears in ANOTHER module's accepted table (run IDs are globally
     unique, and 01b cells are accepted where they are produced);
  3. an accepted run's manifest names it as an upstream;
  4. a git-tracked `_m/combined/` table names it -- those are the citable
     collations, so a run they cite is still being read;
  5. the newest Module 11 figure run, or anything ITS manifest names as
     upstream. This clause is not decoration. Module 01's QC-refresh runs
     (`vmrcatqc-*`) deliberately have no acceptance row of their own -- they
     re-run QC over an already accepted catalog -- yet every Figure 1 panel
     reads them. Without this clause they classify as superseded and a cleanup
     stage would delete the inputs to the manuscript's first figure.

Everything else is `superseded`, `smoke` (named as a smoke or dry run), or
`legacy_v1` (a pre-migration directory, whose retirement MIGRATION_MANIFEST.tsv
governs -- not this file).

`superseded_by` is a HINT, never an instruction: it is the protected run in the
same module sharing the longest name prefix. `review_required` is TRUE wherever
that hint is absent or the run is a diagnostic side-stage whose audit value a
human has to judge. Read those two columns together before deleting anything.

Usage:
    python3 00_shared/deprecated_runs.py [--write]

Without --write it prints a summary and leaves the file alone.
"""

import argparse
import datetime as _dt
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Run IDs look like {stem}-{YYYYMMDD}[-{suffix}], but several legitimate stems
# carry dots (donor-group cells: all_individuals.EA) and the older calibration
# runs carry no date at all. Match generously; the classification never depends
# on parsing a run ID, only the superseded_by hint does.
#
# The suffix is `-?[a-z]`, not `-[a-z]`: this repo mints BOTH `fig-all-20260922-c`
# and `lgv-joint-pve-decision-20260821a`, and 13 live directories use the second
# form. Matching only the hyphenated form truncated those to a shorter ID, so a
# manifest or config citing `...-20260821a` protected `...-20260821` instead and
# left the real run exposed -- under-protection, the one direction that loses
# data. The trailing lookahead stops a match running on into an adjacent word.
RUN_TOKEN = re.compile(
    r"[a-z0-9]+(?:-[A-Za-z0-9_.]+)+-\d{8}(?:-?[a-z])?(?![A-Za-z0-9])")
SMOKE = re.compile(r"smoke|dry|probe|\btest\b", re.I)

# Directories under _m/runs/ that are not runs.
NOT_A_RUN = {"retired", "archive", "logs"}

# Pre-migration trees. AGENTS.md 5.3 governs these through MIGRATION_MANIFEST.tsv,
# so they are recorded here for completeness and explicitly not actionable.
LEGACY_MODULES = {"calibrated-simulation-analysis"}


def sh(*args):
    return subprocess.run(args, capture_output=True, text=True, cwd=REPO).stdout


def accepted_ids(module):
    """Run IDs in a module README's `## Accepted runs` table.

    This reads every table row until the next `##` heading, so a section
    holding a SECOND table collects that table's first column too -- Module 03
    has one keyed by `cell`, which yields tokens like `all_individuals.EA`.
    That is deliberate and must stay. The only use of this set is protection,
    so over-collecting can at worst spare a directory, while under-collecting
    deletes one; a token matching no run directory costs nothing.

    Do not "tighten" this to the first table to match
    `gates.R::read_accepted_runs()`, which is column-validated and correctly
    returns only real rows. The two have opposite failure costs: the gate must
    not admit a bogus acceptance, this must not miss a real one.
    """
    path = os.path.join(REPO, module, "README.md")
    if not os.path.exists(path):
        return set()
    text = open(path, encoding="utf-8", errors="replace").read()
    parts = re.split(r"^##\s+Accepted runs\s*$", text, flags=re.M)
    if len(parts) < 2:
        return set()
    body = re.split(r"^##\s", parts[1], flags=re.M)[0]
    out = set()
    for line in body.splitlines():
        if not line.startswith("|") or set(line) <= set("|- :"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells and cells[0] and cells[0].lower() != "run_id":
            out.add(cells[0])
    return out


def manifest_refs(path):
    if not os.path.exists(path):
        return set()
    return set(RUN_TOKEN.findall(open(path, encoding="utf-8", errors="replace").read()))


def newest_figure_run(modules):
    """The Module 11 run a reader would currently cite: newest SEALED run.

    Sealed matters. An in-flight or abandoned build has no manifest, so it
    protects nothing, and picking it would leave the real figure run's
    upstreams unprotected.
    """
    mod = "11_integrated_manuscript_outputs"
    if mod not in modules:
        return None
    runs = []
    for rid in modules[mod]:
        mf = os.path.join(REPO, mod, "_m", "runs", rid, "manifest.tsv")
        if not os.path.exists(mf):
            continue
        txt = open(mf, encoding="utf-8", errors="replace").read()
        if re.search(r"^finished_at\t\S", txt, flags=re.M):
            runs.append(rid)
    return sorted(runs)[-1] if runs else None


def dir_bytes(path):
    try:
        return int(sh("du", "-sb", path).split()[0])
    except (IndexError, ValueError):
        return 0


def discover_modules():
    """{module: [run_id, ...]} for every module with an `_m/runs/` tree."""
    modules = {}
    for entry in sorted(os.listdir(REPO)):
        runs_dir = os.path.join(REPO, entry, "_m", "runs")
        if not os.path.isdir(runs_dir):
            continue
        modules[entry] = sorted(
            d for d in os.listdir(runs_dir)
            if os.path.isdir(os.path.join(runs_dir, d)) and d not in NOT_A_RUN)
    return modules


def compute_protected(modules):
    """The five protection clauses. Returns (protected_ids, figure_run).

    This is the single source of truth for what may not be deleted, and it is
    deliberately cheap -- no `du`, only reads -- so the cleanup stage can
    recompute it immediately before removing anything rather than trusting a
    DEPRECATED_RUNS.tsv that may have been written weeks earlier. A run that
    became load-bearing after the ledger was generated is invisible to the
    ledger and must not be invisible here.
    """
    accepted = {m: accepted_ids(m) for m in modules}
    protected = set().union(*accepted.values()) if accepted else set()

    # Clauses 3-5.
    for m, ids in accepted.items():
        for rid in ids:
            protected |= manifest_refs(
                os.path.join(REPO, m, "_m", "runs", rid, "manifest.tsv"))
    for f in sh("git", "ls-files", "*/_m/combined/*.tsv").split():
        try:
            protected |= set(RUN_TOKEN.findall(
                open(os.path.join(REPO, f), encoding="utf-8", errors="replace").read()))
        except OSError:
            pass
    # Clause 6: any run named in `config/`. A locked config pins the runs a
    # decision rests on, and those are exactly the runs no accepted result
    # consumes -- so clauses 1-5, which all ask "is this an input to something
    # accepted?", are blind to them by construction.
    #
    # `config/local_genetic_control.yml` is the case that matters:
    # `decision_run_id`, `validation_run_id` and `training_run_id` name the
    # runs behind AGENTS.md 7.2's retirement of `h2_en_calibrated`, and
    # `frozen_model_sha256` pins the model 7.2 requires be applied once. Two of
    # those three classified as `superseded`. They survived the first cleanup
    # tranche only because no successor could be guessed for them, which is a
    # naming heuristic, not knowledge that the project depends on them.
    for root, _, files in os.walk(os.path.join(REPO, "config")):
        for name in files:
            try:
                protected |= set(RUN_TOKEN.findall(
                    open(os.path.join(root, name), encoding="utf-8",
                         errors="replace").read()))
            except OSError:
                pass

    fig = newest_figure_run(modules)
    if fig:
        protected.add(fig)
        protected |= manifest_refs(os.path.join(
            REPO, "11_integrated_manuscript_outputs", "_m", "runs", fig, "manifest.tsv"))
    return protected, fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="write DEPRECATED_RUNS.tsv (otherwise summarise only)")
    args = ap.parse_args()

    modules = discover_modules()
    protected, fig = compute_protected(modules)

    def hint(module, rid):
        """Longest-prefix protected run in the same module."""
        best, best_len = "", 0
        for cand in modules[module]:
            if cand == rid or cand not in protected:
                continue
            n = len(os.path.commonprefix([cand, rid]))
            if n > best_len:
                best, best_len = cand, n
        # A handful of shared characters is not evidence of succession.
        return best if best_len >= max(8, len(rid) // 2) else ""

    rows, skipped = [], 0
    for m in modules:
        for rid in modules[m]:
            if rid in protected:
                skipped += 1
                continue
            if m in LEGACY_MODULES:
                status, reason, review = ("legacy_v1",
                                          "pre-migration tree; retirement governed by "
                                          "MIGRATION_MANIFEST.tsv, not this file", "TRUE")
            elif SMOKE.search(rid):
                status, reason, review = ("smoke",
                                          "smoke or dry run; never intended for acceptance",
                                          "FALSE")
            else:
                status, reason, review = ("superseded",
                                          "no accepted result depends on this run", "FALSE")
            sb = hint(m, rid)
            if status == "superseded" and not sb:
                review = "TRUE"
                reason += "; no successor identified"
            path = f"{m}/_m/runs/{rid}"
            rows.append({
                "module": m, "run_id": rid, "path": path, "status": status,
                "bytes": str(dir_bytes(os.path.join(REPO, path))),
                "superseded_by": sb, "reason": reason,
                "review_required": review,
                "deprecated_on": _dt.date.today().isoformat(),
                "cleanup_status": "pending",
            })

    cols = ["module", "run_id", "path", "status", "bytes", "superseded_by",
            "reason", "review_required", "deprecated_on", "cleanup_status"]

    # Carry forward what the cleanup stage already did. This survey only sees
    # directories that still exist, so a deleted run would silently vanish from
    # the ledger on the next --write and take the record of its deletion with
    # it -- leaving no answer to "what was removed, and when". Rows already
    # marked `deleted:` are preserved verbatim, and a previously recorded
    # `cleanup_status` is kept for any run still on disk so a human decision
    # (`keep`, `hold`, ...) is not reset by re-running the survey.
    previous = {}
    out = os.path.join(REPO, "DEPRECATED_RUNS.tsv")
    if os.path.exists(out):
        with open(out, encoding="utf-8") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if line.strip():
                    r = dict(zip(head, line.rstrip("\n").split("\t")))
                    previous[r.get("run_id", "")] = r

    live = {r["run_id"] for r in rows}
    for r in rows:
        prior = previous.get(r["run_id"], {}).get("cleanup_status", "")
        if prior and prior != "pending":
            r["cleanup_status"] = prior
    for rid, r in previous.items():
        if rid not in live and r.get("cleanup_status", "").startswith("deleted"):
            rows.append({c: r.get(c, "") for c in cols})

    rows.sort(key=lambda r: (r["module"], r["run_id"]))

    # Summarise only what is still on disk. The carried-forward `deleted:` rows
    # are history, and counting their bytes here would report space that has
    # already been reclaimed as though it were still waiting to be.
    pending = [r for r in rows if not r["cleanup_status"].startswith("deleted")]
    gone = [r for r in rows if r["cleanup_status"].startswith("deleted")]

    total = sum(int(r["bytes"]) for r in pending)
    by = {}
    for r in pending:
        s = by.setdefault(r["status"], [0, 0])
        s[0] += 1
        s[1] += int(r["bytes"])
    print(f"protected (kept, not listed): {skipped}")
    for s, (n, b) in sorted(by.items()):
        print(f"  {s:<12} {n:4d} dirs  {b / 1e9:7.2f} GB")
    print(f"  {'TOTAL':<12} {len(pending):4d} dirs  {total / 1e9:7.2f} GB")
    print(f"  review_required: {sum(1 for r in pending if r['review_required'] == 'TRUE')}")
    if gone:
        freed = sum(int(r["bytes"]) for r in gone)
        print(f"  already deleted: {len(gone):4d} dirs  {freed / 1e9:7.2f} GB "
              f"(history, retained in this file)")
    if fig:
        print(f"current figure run protected: {fig}")

    if not args.write:
        print("\n(dry run; pass --write to update DEPRECATED_RUNS.tsv)")
        return 0

    tmp = out + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    os.replace(tmp, out)
    print(f"\nwrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
