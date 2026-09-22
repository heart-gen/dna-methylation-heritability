#!/usr/bin/env python3
"""Delete run directories that `deprecated_runs.py` marked and a human cleared.

This is the removal half of the two-stage retirement AGENTS.md 3 requires. It is
a separate script from the marking stage on purpose: marking is a survey that
can be re-run freely, deletion is irreversible, and the two must not share a
command line where a stray flag turns one into the other.

Irreversible is meant literally here. Run directories are NOT tracked in git --
`git ls-files` returns nothing for any of the 137 marked paths -- so there is no
commit, stash or reflog behind them. AGENTS.md 3 allows historical recovery from
"a tagged Git commit, archived release, or existing Zenodo record"; for these
trees none of those exists, and the only recovery path is re-running the
analysis from `_h/` plus locked configuration, which AGENTS.md 5.2 guarantees is
possible. Deleting a run therefore costs the compute to rebuild it, not the
ability to rebuild it -- but it does cost that, so the guards below are not
ceremony.

What it will delete: rows whose `review_required` is FALSE and whose
`cleanup_status` is still `pending`. Those are the runs with an identified
successor and the smoke runs. Everything a human still has to judge --
`review_required = TRUE`, which includes every `legacy_v1` row governed by
MIGRATION_MANIFEST.tsv -- is left alone and needs a deliberate flip of that
column to become eligible.

Four guards run before anything is removed:

  1. protection is RECOMPUTED live via `deprecated_runs.compute_protected()`,
     not read from the ledger. A ledger is a snapshot; a run that became
     load-bearing after it was written -- a new acceptance row, a new figure
     run -- is invisible to the snapshot and must not be invisible here. Any
     overlap is a hard error that deletes nothing, because it means the ledger
     and the repository disagree and a human should find out why;
  2. each path must resolve inside this repository, sit under an `_m/runs/`
     tree, and end in exactly its own `run_id`. A corrupted or hand-edited row
     cannot steer `rm -rf` somewhere else;
  3. the recorded `bytes` is compared against the directory's current size, and
     a mismatch beyond a small tolerance stops that row. A run that changed
     size since it was surveyed is a run something still writes to;
  4. provenance is archived first (see below), and a row whose archive step
     fails is not deleted.

A sealed run survives this script, by design. `close_run()` leaves a finished
run's directory without its write bit (AGENTS.md 5.2 immutability), and a
directory with no write bit will not let its entries be unlinked, so `rmtree`
reports "Permission denied" and the row stays `pending`. That is the seal
working, not a bug, and this script deliberately does NOT clear it: forcing a
recursive delete through a permission the repository set on purpose is a
decision a person should make per run, not a default buried in a cleanup tool.
Such rows are reported as failures so they stay visible. To remove one, clear
its seal explicitly first (`chmod -R u+w <path>`) and re-run.

Provenance outlives the data. Before a directory is removed, its `manifest.tsv`
and `output_checksums.tsv` are copied into `_deleted_run_provenance/` at the
repository root. Those two files are what AGENTS.md 9 actually requires a run to
carry -- git commit, config checksums, upstream run IDs, output checksums -- so
keeping them means a deleted run can still be audited, and anything rebuilt in
its place can be checksum-verified against what it replaced.

The archive is not free but it is cheap: the first 72 deletions kept 43 MB
against 46 GB removed, roughly 0.1%. Almost all of that is three Module 09
`output_checksums.tsv` files of 11-16 MB each, one line per output in runs that
emitted hundreds of thousands. If the ratio ever stops looking cheap, drop
`output_checksums.tsv` from KEEP before `manifest.tsv` -- checksums can be
recomputed from a rebuild, but a run's identity cannot be reconstructed at all.

Usage:
    python3 00_shared/cleanup_deprecated_runs.py            # plan only
    python3 00_shared/cleanup_deprecated_runs.py --execute  # delete

Without --execute it prints the plan and removes nothing.
"""

import argparse
import datetime as _dt
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deprecated_runs as dr  # noqa: E402

REPO = dr.REPO
LEDGER = os.path.join(REPO, "DEPRECATED_RUNS.tsv")
PROVENANCE = os.path.join(REPO, "_deleted_run_provenance")

# Files worth keeping after the run itself is gone.
KEEP = ("manifest.tsv", "output_checksums.tsv")

# `du` and the recorded size can differ slightly through sparse files and
# filesystem accounting. A real write is far larger than this.
SIZE_TOLERANCE = 0.02


def read_ledger():
    with open(LEDGER, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return header, [dict(zip(header, line.rstrip("\n").split("\t")))
                        for line in fh if line.strip()]


def safe_path(row):
    """Resolve a ledger path, refusing anything that is not this run's dir."""
    rid = row["run_id"]
    if not rid or "/" in rid or ".." in rid:
        return None, f"implausible run_id {rid!r}"
    full = os.path.realpath(os.path.join(REPO, row["path"]))
    if not full.startswith(os.path.realpath(REPO) + os.sep):
        return None, "resolves outside the repository"
    if f"{os.sep}_m{os.sep}runs{os.sep}" not in full + os.sep:
        return None, "not under an _m/runs/ tree"
    if os.path.basename(full) != rid:
        return None, "basename does not match run_id"
    return full, None


def archive_provenance(full, row):
    """Copy the run's provenance files out before the tree is removed."""
    dest = os.path.join(PROVENANCE, row["module"], row["run_id"])
    saved = []
    for name in KEEP:
        src = os.path.join(full, name)
        if not os.path.exists(src):
            continue
        os.makedirs(dest, exist_ok=True)
        shutil.copy2(src, os.path.join(dest, name))
        saved.append(name)
    return saved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true",
                    help="actually delete (default: print the plan only)")
    args = ap.parse_args()

    header, rows = read_ledger()
    if "cleanup_status" not in header:
        print("ERROR: ledger has no cleanup_status column", file=sys.stderr)
        return 1

    modules = dr.discover_modules()
    protected, fig = dr.compute_protected(modules)
    print(f"protection recomputed: {len(protected)} run IDs protected")
    print(f"current figure run:    {fig}\n")

    eligible = [r for r in rows
                if r["review_required"] == "FALSE"
                and r["cleanup_status"] == "pending"]
    held = len(rows) - len(eligible)

    # Guard 1. A ledger/repository disagreement stops the whole stage.
    clash = sorted(r["run_id"] for r in eligible if r["run_id"] in protected)
    if clash:
        print("ERROR: these runs are marked deletable but are PROTECTED now:",
              file=sys.stderr)
        for rid in clash:
            print(f"  {rid}", file=sys.stderr)
        print("\nThe ledger is stale. Re-run deprecated_runs.py --write and "
              "re-read it before deleting anything.", file=sys.stderr)
        return 1

    plan, refused = [], []
    for row in eligible:
        full, why = safe_path(row)
        if why:
            refused.append((row["run_id"], why))
            continue
        if not os.path.isdir(full):
            refused.append((row["run_id"], "directory is already gone"))
            continue
        # Guard 3.
        recorded, actual = int(row["bytes"]), dr.dir_bytes(full)
        if recorded > 0 and abs(actual - recorded) / recorded > SIZE_TOLERANCE:
            refused.append((row["run_id"],
                            f"size changed since survey: {recorded} -> {actual}"))
            continue
        plan.append((row, full, actual))

    total = sum(b for _, _, b in plan)
    print(f"eligible to delete : {len(plan):4d} runs  {total / 1e9:7.2f} GB")
    print(f"refused by guards  : {len(refused):4d}")
    print(f"held for review    : {held:4d} runs "
          f"(review_required=TRUE or already handled)\n")
    for rid, why in refused:
        print(f"  REFUSED {rid}: {why}")
    if refused:
        print()

    by_mod = {}
    for row, _, b in plan:
        m = by_mod.setdefault(row["module"], [0, 0])
        m[0] += 1
        m[1] += b
    for m, (n, b) in sorted(by_mod.items()):
        print(f"  {m:<38} {n:3d} runs  {b / 1e9:7.2f} GB")

    if not args.execute:
        print("\n(plan only; pass --execute to delete)")
        return 0

    print(f"\nDeleting {len(plan)} run directories...")
    today = _dt.date.today().isoformat()
    deleted, failed = 0, []
    for row, full, _ in plan:
        try:
            saved = archive_provenance(full, row)  # guard 4
        except OSError as exc:
            failed.append((row["run_id"], f"provenance archive failed: {exc}"))
            continue
        try:
            shutil.rmtree(full)
        except OSError as exc:
            failed.append((row["run_id"], f"rmtree failed: {exc}"))
            continue
        row["cleanup_status"] = f"deleted:{today}"
        deleted += 1
        print(f"  deleted {row['run_id']} (kept {', '.join(saved) or 'nothing'})")

    for rid, why in failed:
        print(f"  FAILED {rid}: {why}", file=sys.stderr)

    tmp = LEDGER + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in header) + "\n")
    os.replace(tmp, LEDGER)

    print(f"\ndeleted {deleted} runs; ledger updated")
    print(f"provenance kept in {os.path.relpath(PROVENANCE, REPO)}/")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
