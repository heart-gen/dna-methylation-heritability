#!/usr/bin/env python
"""Restore the annotation name on an LD score file's score column.

This LDSC build writes the LD-score column of a single ``--thin-annot`` file as
a bare ``L2``, discarding the annotation name carried in the ``.annot.gz``
header. That name is what ``ldsc.py --h2`` turns into the ``Category`` label of
the ``.results`` file, so without this rewrite the annotation arrives
downstream as the positional label ``L2_1`` and stage 06 cannot identify it by
name. Restoring the name keeps every ``.results`` file self-describing and
keeps no stage dependent on row order.

Idempotent: re-running on an already-labelled file is a no-op.
"""
import argparse
import gzip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot", required=True)
    ap.add_argument("--ldscore", required=True)
    args = ap.parse_args()

    with gzip.open(args.annot, "rt") as fh:
        names = fh.readline().rstrip("\n").split("\t")
    if len(names) != 1:
        raise SystemExit(
            "expected a single --thin-annot column in %s, found %d: %s"
            % (args.annot, len(names), names)
        )
    want = names[0] + "L2"

    with gzip.open(args.ldscore, "rt") as fh:
        lines = fh.readlines()
    header = lines[0].rstrip("\n").split("\t")
    if header[-1] == want:
        return
    if header[-1] != "L2":
        raise SystemExit(
            "unexpected LD score column name %r in %s, refusing to rewrite"
            % (header[-1], args.ldscore)
        )
    header[-1] = want
    lines[0] = "\t".join(header) + "\n"
    with gzip.open(args.ldscore, "wt") as fh:
        fh.writelines(lines)


if __name__ == "__main__":
    main()
