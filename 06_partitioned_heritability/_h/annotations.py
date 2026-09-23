"""06_partitioned_heritability -- which annotations enter the S-LDSC model.

Two annotations, in this column order:

    1. VMR_TESTED                 binary membership in the tested VMR universe
    2. LOCAL_SNP_CONTRIBUTION_Z   the continuous within-cell relative score

**Why two and not one.** `local_snp_contribution_score_z` is a standardized
within-cell midrank percentile, so it is symmetric about zero *by construction*
and zero is the value of a median-ranked VMR. A single-column annotation
therefore hands a SNP in no VMR and a SNP in a median-ranked VMR the same
number, and with no membership term in the model there is nothing to absorb a
VMR-versus-genome difference. tau then blurs the within-VMR gradient the module
asks about with a membership effect it never meant to test, which is why the
one-column result was not merely negative but uninterpretable.

With `VMR_TESTED` in the model, tau on the score is conditional on membership
and estimates the within-VMR gradient. tau on `VMR_TESTED` is a second,
separately interesting quantity -- are VMRs enriched at all -- reported
descriptively and outside the frozen FDR family, so that family's size and its
already-computed q-values do not change.

See `writing-notes/ISSUE_06_two_annotation_model.md`.

**The order is part of the contract.** `.annot.gz` columns, the score columns of
`.l2.ldscore.gz`, and the entries of `.l2.M` / `.l2.M_5_50` are all positional,
so every stage has to agree on it. Downstream identification is nonetheless by
NAME and never by position: LDSC labels each LD-score column `<NAME>L2` when
there is more than one annotation, `ldscore_fromlist` appends `_<file index>`,
and `Category` in the `.results` file inherits that label. Membership is placed
first so the primary annotation remains the last column, matching the
baselineLD-then-ours layout the rest of the module describes.

Neither name may be a prefix of the other: stage 06 identifies rows by name
prefix, and a nested pair would make that match ambiguous. That is asserted at
import, not merely documented.
"""
from __future__ import annotations

import re

MEMBERSHIP_ANNOT = "VMR_TESTED"
SCORE_ANNOT = "LOCAL_SNP_CONTRIBUTION_Z"

#: Column order of the `.annot.gz` file. Positional; see the module docstring.
ANNOT_COLUMNS = (MEMBERSHIP_ANNOT, SCORE_ANNOT)

#: The annotation carrying the module's primary hypothesis.
PRIMARY_ANNOT = SCORE_ANNOT

#: Role per annotation. Downstream R stages key on the role rather than on the
#: name, so the name lives in exactly one place (here).
ANNOT_ROLE = {MEMBERSHIP_ANNOT: "membership", SCORE_ANNOT: "score"}

#: Enrichment is a ratio of heritability share to annotation share. For a signed
#: continuous annotation the denominator is a signed sum, so the ratio is not
#: interpretable; for a binary indicator it is the usual S-LDSC enrichment.
ENRICHMENT_INTERPRETABLE = {MEMBERSHIP_ANNOT: True, SCORE_ANNOT: False}

#: `<NAME>L2` optionally suffixed `_<index>` by ldsc's sideways concatenation.
_LDSCORE_SUFFIX = re.compile(r"^L2(_\d+)?$")

#: The config key that pins the source column of the continuous annotation. The
#: annotation's own label is fixed here instead of derived from it, so accepted
#: runs' `.results` category names stay comparable across runs.
SCORE_SOURCE_COLUMN = "local_snp_contribution_score_z"


def _assert_names_unambiguous() -> None:
    for a in ANNOT_COLUMNS:
        for b in ANNOT_COLUMNS:
            if a is not b and b.startswith(a):
                raise SystemExit(
                    f"Annotation names {a!r} and {b!r} are nested; stage 06 "
                    "identifies .results rows by name prefix and could not "
                    "tell them apart.")


_assert_names_unambiguous()


def ldscore_columns() -> list:
    """Names LDSC gives the LD-score columns of a multi-annotation file."""
    return [name + "L2" for name in ANNOT_COLUMNS]


def check_annot_header(names, source: str) -> None:
    """Fail closed unless `names` is exactly the declared contract, in order.

    A single-column header is one of the things this refuses: it is the
    one-annotation model the two-annotation construction replaced, and letting
    it through would silently produce an uninterpretable tau.
    """
    got = list(names)
    if got != list(ANNOT_COLUMNS):
        raise SystemExit(
            f"{source}: annotation columns are {got}, expected exactly "
            f"{list(ANNOT_COLUMNS)} in that order. Module 06 runs a "
            "two-annotation model (membership + continuous score); a "
            "one-annotation file is the construction this replaced and is "
            "refused rather than analysed.")


def assert_score_source_column(score_column: str) -> None:
    if score_column != SCORE_SOURCE_COLUMN:
        raise SystemExit(
            f"config annotation.score_column is {score_column!r}; Module 06 "
            f"requires {SCORE_SOURCE_COLUMN!r} (AGENTS.md 3).")


def find_results_row(df, name: str):
    """Return the single `.results` row for annotation `name`.

    Matching is by name plus an exact LD-score suffix, never by position: a
    baselineLD version change must not be able to shift which row is read, and
    a loose prefix match must not be able to pick up a neighbouring category.
    """
    categories = df.iloc[:, 0].astype(str)
    keep = [c.startswith(name) and bool(_LDSCORE_SUFFIX.match(c[len(name):]))
            for c in categories]
    hit = df[keep]
    if hit.empty:
        raise SystemExit(
            f"no row named {name}L2* in the .results file. "
            f"Found: {list(categories)[-4:]}")
    if len(hit) > 1:
        raise SystemExit(
            f"{len(hit)} rows match {name}L2*: {list(hit.iloc[:, 0])}")
    return hit.iloc[0]
