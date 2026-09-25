#!/usr/bin/env python
"""Check -- and if needed restore -- the annotation names on an LD score file.

``ldsc.py --l2`` names the LD-score columns ``<ANNOT>L2`` when the annotation
file has more than one column, but writes a bare ``L2`` when it has exactly one
(``ldsc.py``: ``if n_annot == 1: ldscore_colnames = [col_prefix+scale_suffix]``).
That name is what ``ldsc.py --h2`` turns into the ``Category`` label of the
``.results`` file, so under the old one-annotation model the annotation arrived
downstream as the positional label ``L2_1`` and stage 06 could not identify it
by name. This stage restored the name.

Module 06 now runs a two-annotation model, for which LDSC labels the columns
itself, so on the normal path this stage has nothing to rewrite. It is kept, and
kept in the chain, for two reasons:

* it is the point where the annotation contract is checked against what actually
  reached the LD scores, before eight S-LDSC regressions are launched off it;
* the single-column rewrite remains implemented, so the stage tells a bare
  ``L2`` apart from an unexpected name rather than conflating them.

Fail-closed on both sides. The ``.annot.gz`` header must be exactly the declared
annotation set, in the declared order -- a one-column file is refused here
rather than silently analysed as the one-annotation model this replaced -- and
the LD-score file's trailing columns must be either already-correct names or the
bare ``L2`` of the single-annotation case. Anything else aborts.

Idempotent: re-running on an already-labelled file is a no-op.
"""
import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotations import check_annot_header, ldscore_columns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", required=True)
    ap.add_argument("--ldscore", required=True)
    args = ap.parse_args()

    with gzip.open(args.annot, "rt") as fh:
        names = fh.readline().rstrip("\n").split("\t")
    check_annot_header(names, args.annot)
    want = ldscore_columns()
    n = len(want)

    with gzip.open(args.ldscore, "rt") as fh:
        lines = fh.readlines()
    header = lines[0].rstrip("\n").split("\t")
    if len(header) < n:
        raise SystemExit(
            "%s has %d columns, too few for %d annotations plus CHR/SNP/BP"
            % (args.ldscore, len(header), n))
    tail = header[-n:]
    if tail == want:
        return
    # The only other shape this stage accepts: the single-annotation case, where
    # LDSC drops the name entirely. Naming it by position is safe only because
    # there is exactly one column to name.
    if n == 1 and tail == ["L2"]:
        header[-1] = want[0]
        lines[0] = "\t".join(header) + "\n"
        with gzip.open(args.ldscore, "wt") as fh:
            fh.writelines(lines)
        return
    raise SystemExit(
        "unexpected LD score column name(s) %r in %s (expected %r), refusing "
        "to rewrite" % (tail, args.ldscore, want))


if __name__ == "__main__":
    main()
