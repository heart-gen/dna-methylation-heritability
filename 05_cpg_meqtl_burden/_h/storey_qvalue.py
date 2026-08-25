#!/usr/bin/env python
"""Storey q-values for 02b_combine_meqtl.R, via py_qvalue.

Runs under the `genomics` environment, which is where py_qvalue is installed;
the R stage lives in `epigenomics`, which has neither the Bioconductor `qvalue`
package nor py_qvalue. Rather than add a package to a shared environment, the
R stage shells out to this script.

`lfdr_out=False` is not an optimisation, it is required: py_qvalue's local-FDR
branch does not return in any reasonable time even on a few thousand p-values,
and the module needs q-values only.

Reads one p-value per line on stdin (no header), writes `qvalue` on the first
line as `pi0<TAB>value` and then one q-value per line, in input order.
"""
import sys

import numpy as np
import py_qvalue


def main() -> int:
    p = np.loadtxt(sys.stdin, dtype=float, ndmin=1)
    if p.size == 0:
        print("ERROR\tno p-values on stdin", file=sys.stderr)
        return 1
    if not np.all(np.isfinite(p)) or p.min() < 0 or p.max() > 1:
        print("ERROR\tp-values must be finite and in [0, 1]", file=sys.stderr)
        return 1

    res = py_qvalue.qvalue(p, lfdr_out=False)
    q = np.asarray(res["qvalues"], dtype=float)
    if q.shape != p.shape:
        print(f"ERROR\tgot {q.shape} q-values for {p.shape} p-values",
              file=sys.stderr)
        return 1

    # Storey q-values are monotone in p and bounded above by pi0. Checking here
    # means the R stage never has to trust an unvalidated 0.1.0 dependency.
    order = np.argsort(p, kind="stable")
    pi0 = float(res["pi0"])
    if np.any(np.diff(q[order]) < -1e-9):
        print("ERROR\tq-values are not monotone in p", file=sys.stderr)
        return 1
    if not (0 < pi0 <= 1) or q.max() > pi0 + 1e-9 or q.min() < 0:
        print(f"ERROR\timplausible pi0={pi0} or q range "
              f"[{q.min()}, {q.max()}]", file=sys.stderr)
        return 1

    out = sys.stdout
    out.write(f"pi0\t{pi0!r}\n")
    np.savetxt(out, q, fmt="%.17g")
    return 0


if __name__ == "__main__":
    sys.exit(main())
