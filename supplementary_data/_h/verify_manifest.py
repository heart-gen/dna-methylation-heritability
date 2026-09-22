#!/usr/bin/env python3
"""Check supplementary_data_manifest.tsv against what is actually on disk.

The manifest declares WHERE each supplementary-data file lives. Nothing is
copied or symlinked into `supplementary_data/`, so the declared locations are
the only thing keeping the list honest -- this script is what verifies them.

  python3 supplementary_data/_h/verify_manifest.py           # report only
  python3 supplementary_data/_h/verify_manifest.py --write    # refresh sizes

Exits non-zero when a location that is supposed to exist resolves to nothing,
so it is safe to run before staging a Zenodo upload.
"""
import csv, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MANIFEST = os.path.join(ROOT, "supplementary_data", "supplementary_data_manifest.tsv")

## A location may be written with one brace group, `dir/{a,b,c}.tsv`, because
## listing nine sibling tables on nine rows would bury the one thing a reader
## wants from the row. Expansion happens here, before globbing.
def expand_braces(pat):
    m = re.search(r"\{([^{}]*)\}", pat)
    if not m:
        return [pat]
    pre, post = pat[:m.start()], pat[m.end():]
    out = []
    for alt in m.group(1).split(","):
        out.extend(expand_braces(pre + alt + post))
    return out

def resolve(location):
    """Return the list of real files a location names, relative to ROOT."""
    import glob as g
    hits = []
    for pat in expand_braces(location):
        for p in g.glob(os.path.join(ROOT, pat)):
            if os.path.isdir(p):
                for dirpath, _, names in os.walk(p):
                    hits += [os.path.join(dirpath, n) for n in names]
            elif os.path.isfile(p):
                hits.append(p)
    return sorted(set(hits))

def tracked(paths):
    """How many of these paths git tracks. `_m/runs/` is gitignored by design,
    so a run output answering 0 here is correct, not a problem."""
    if not paths:
        return 0
    rel = [os.path.relpath(p, ROOT) for p in paths]
    n = 0
    for i in range(0, len(rel), 400):
        r = subprocess.run(["git", "-C", ROOT, "ls-files", "-z", "--"] + rel[i:i + 400],
                           capture_output=True, text=True)
        n += len([x for x in r.stdout.split("\0") if x])
    return n

def human(b):
    for unit in ("B", "K", "M", "G"):
        if b < 1024 or unit == "G":
            return f"{b:.0f}{unit}" if unit == "B" else f"{b:.1f}{unit}"
        b /= 1024

def main():
    write = "--write" in sys.argv
    with open(MANIFEST) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    hdr = list(rows[0].keys())

    problems, by_sd = [], {}
    for row in rows:
        paths = resolve(row["location"])
        n, b = len(paths), sum(os.path.getsize(p) for p in paths)
        row["n_files"], row["bytes"] = str(n), str(b)
        row["in_git"] = "yes" if n and tracked(paths) == n else \
                        "partial" if tracked(paths) else "no"
        ## A row that claims to be deposit-ready and resolves to nothing is the
        ## failure this script exists to catch. `pending_build` has no files yet
        ## on purpose, and an excluded directory may have been cleaned off Quest.
        if n == 0 and row["status"] not in ("pending_build", "excluded"):
            problems.append(f"  {row['location']}  [{row['status']}]")
        ## Statuses that contribute nothing to a Zenodo upload: retired or
        ## blocked outputs, external data, Quest-only bulk, a build that has not
        ## happened, and `partial_deposit`, where only a named subset ships and
        ## the resolved size would overstate it by two orders of magnitude.
        if row["status"] in ("ready", "needs_deid", "pending_acceptance",
                             "transform_on_deposit"):
            k = (row["sd"], row["title"])
            by_sd[k] = by_sd.get(k, 0) + b

    print(f"{'SD':>3}  {'bytes':>8}  title")
    for (sd, title), b in sorted(by_sd.items(), key=lambda x: int(x[0][0])):
        print(f"{sd:>3}  {human(b):>8}  {title}")
    print(f"\ndeposit total: {human(sum(by_sd.values()))}")
    print(f"rows: {len(rows)}   locations resolving to 0 files: {len(problems)}")
    if problems:
        print("\nunresolved:")
        print("\n".join(problems))

    if write:
        with open(MANIFEST, "w", newline="") as fh:
            w = csv.DictWriter(fh, hdr, delimiter="\t", lineterminator="\n",
                               quoting=csv.QUOTE_NONE, escapechar=None)
            w.writeheader()
            w.writerows(rows)
        print(f"\nrewrote {os.path.relpath(MANIFEST, ROOT)}")
    return 1 if problems else 0

if __name__ == "__main__":
    sys.exit(main())
