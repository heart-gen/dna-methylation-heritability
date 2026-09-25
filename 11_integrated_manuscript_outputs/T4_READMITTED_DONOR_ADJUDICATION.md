# T4 — the four flagged readmitted donors: adjudication packet

Prepared 2026-09-25 against the accepted run `fig-all-20260922-c`. Every number
below is read from
`11_integrated_manuscript_outputs/_m/runs/fig-all-20260922-c/tables/qc_cross_region_concordance.tsv`
(320 donor x region-pair rows, 120 donors) and
`.../tables/qc_swap_candidates.tsv`. Nothing is carried over from TASKS.md.

**DECIDED 2026-09-25 by the PI: all eight donors are cleared and
`config/cohorts.yml:300` is affirmed.** The decision and its grounds are in
section 6; the affirmation is recorded in the `sample_blacklist` comment block of
`config/cohorts.yml`. T4 is closed. No donor is excluded, so no module reruns and
`vmr_set_id` is unchanged in every region.

This document is the record of what was decided and why. It is retained rather
than deleted because the reasoning in section 4 is the substance of the decision,
and because sections 3 and 4 supersede numbers that are still quoted elsewhere.

**Why it is urgent rather than merely open.** All eight donors are in the AA
analysis set (they appear in `01_vmr_catalog/_m/runs/vmrcat-AA-*-20260816/plink_format/*/*.fam`).
Excluding any one changes the donor set, which re-derives Module 01, which changes
`vmr_set_id`, which invalidates **every** downstream module including the three
Module 05 runs now awaiting acceptance. This is the only open item in the project
with that reach, so it should be settled **before** the rerun cascade, not after.

---

## 1. What the screen is, in its own words

From `11_integrated_manuscript_outputs/_h/05_qc_sample_integrity.R:12-22`:

- The test is `z_self < 0`: the donor's methylation profile resembles the average
  *other* donor more than it resembles its own other region.
- **Module 01 masks C->T SNP CpGs** (`00_prepare.R remove_ct_snps`) — exactly the
  genotype-driven sites a methylation fingerprint depends on. What is left is
  "a real but weak identity signal": 45-65% of donors are their own single best
  match out of ~150 candidates, against a ~0.7% chance rate.
- The script's own conclusion: *"far above chance but not clean enough to
  adjudicate an individual donor, so this is a SCREEN that nominates candidates
  for follow-up, not a verdict."*
- It is deliberately **not** a genotype-concordance check: genotypes come from one
  draw per donor and are shared across regions by construction, so they cannot
  detect a swap introduced at the methylation assay.
- One autosome (chr22), up to 5,000 CpGs. Nothing downstream consumes it.

So the screen cannot convict a donor on its own. It can, however, be read for
*internal coherence*, which is what section 4 does and which turns out to settle
the question.

## 2. The eight readmitted donors, every region-pair row

Source: `config/cohorts.yml:308-310` — `Br1249 Br1303 Br1371 Br1552 Br1693
Br1700 Br1883` from the DLPFC list and `Br1927` from the hippocampus list.
`rank_self` is the donor's own rank among `n_candidates` possible matches; rank 1
means the donor's best match in the whole cohort is itself.

| donor | pair | rank_self | n_cand | z_self | flagged |
|---|---|---|---|---|---|
| Br1249 | caud-dlpf | 4 | 103 | +1.567 | |
| Br1249 | caud-hipp | **1** | 102 | **+4.820** | |
| Br1249 | dlpf-hipp | 75 | 115 | −0.521 | **TRUE** |
| Br1303 | caud-dlpf | 22 | 103 | +0.767 | |
| Br1303 | caud-hipp | 2 | 102 | +1.943 | |
| Br1303 | dlpf-hipp | 8 | 115 | +1.593 | |
| Br1371 | caud-dlpf | 16 | 103 | +1.130 | |
| Br1371 | caud-hipp | **1** | 102 | +3.478 | |
| Br1371 | dlpf-hipp | 53 | 115 | +0.175 | |
| Br1552 | caud-dlpf | 3 | 103 | +1.674 | |
| Br1552 | caud-hipp | **1** | 102 | **+4.924** | |
| Br1552 | dlpf-hipp | 45 | 115 | +0.345 | |
| Br1693 | caud-dlpf | 49 | 103 | −0.187 | **TRUE** |
| Br1693 | caud-hipp | 22 | 102 | +0.848 | |
| Br1693 | dlpf-hipp | 6 | 115 | +1.916 | |
| Br1700 | caud-dlpf | 22 | 103 | +0.730 | |
| Br1700 | caud-hipp | **1** | 102 | +3.525 | |
| Br1700 | dlpf-hipp | 87 | 115 | −0.770 | **TRUE** |
| Br1883 | caud-dlpf | 19 | 103 | +0.908 | |
| Br1883 | caud-hipp | **1** | 102 | +4.388 | |
| Br1883 | dlpf-hipp | 42 | 115 | +0.394 | |
| Br1927 | caud-dlpf | **1** | 103 | **+4.891** | |
| Br1927 | caud-hipp | 2 | 102 | +1.896 | |
| Br1927 | dlpf-hipp | 54 | 115 | +0.094 | |

**Three of eight are flagged, not four.** `Br1927` — one of the four T4 names — is
`swap_candidate = FALSE` in all three of its pairs and is its own best match in
caudate-DLPFC out of 103 candidates (z = +4.891). It should be dropped from the
question entirely.

| donor | flagged in | best rank achieved | verdict |
|---|---|---|---|
| Br1249 | 1 of 3 pairs | 1 | flagged |
| Br1693 | 1 of 3 pairs | 6 | flagged |
| Br1700 | 1 of 3 pairs | 1 | flagged |
| Br1303 | 0 of 3 | 2 | clean |
| Br1371 | 0 of 3 | 1 | clean |
| Br1552 | 0 of 3 | 1 | clean |
| Br1883 | 0 of 3 | 1 | clean |
| Br1927 | 0 of 3 | 1 | clean |

## 3. The association, both units, recomputed

T4 records p = 0.0064 / OR 9.8. That is reproduced exactly by a **4-of-8 vs
10-of-112** table — 14 flagged donors of 120. The accepted run flags **9 of 120
(7.5%)** and **3 of the 8**. T4's table does not describe this run; the run it was
written against (`fig-all-20260826`, `-a`) was deleted 2026-09-22, so the two
screens cannot be diffed.

**Unit = donor** (a donor is flagged if any of its pairs is):

|  | flagged | clean |
|---|---|---|
| readmitted | 3 | 5 |
| other | 6 | 106 |

Fisher two-sided **p = 0.0137**, OR 10.16, 95% CI **1.28 – 69.48**.
Flag rate 37.5% readmitted vs 5.4% other.

**Unit = donor x region-pair:**

|  | flagged | clean |
|---|---|---|
| readmitted | 3 | 21 |
| other | 8 | 288 |

Fisher two-sided **p = 0.0411**, OR 5.09, 95% CI **0.81 – 23.34** — the interval
**includes 1**.

Neither unit was prespecified. **Do not quote a p without its unit.** On three
events the confidence intervals span one to nearly two orders of magnitude, and
the two units disagree about whether the interval excludes the null.

## 4. The decisive test: is any flag consistent with an actual swap?

A mislabelled tube is a property of **one region's sample**. If donor X's DLPFC
tube is really donor Y's, then **both** pairs involving DLPFC must break —
caudate-DLPFC and DLPFC-hippocampus — while caudate-hippocampus stays clean. This
is a strong, falsifiable prediction and the table can test it donor by donor.

| donor | flagged pairs | best single-tube explanation | pairs that break under it | coherent? | readmitted |
|---|---|---|---|---|---|
| Br1134 | dlpf-hipp, caud-hipp | hippocampus | **2 of 2** | **YES** | no |
| Br1459 | caud-dlpf, dlpf-hipp | DLPFC | **2 of 2** | **YES** | no |
| Br1615 | dlpf-hipp | either | 1 of 2 | no | no |
| Br1185 | dlpf-hipp | either | 1 of 2 | no | no |
| Br1054 | caud-hipp | either | 1 of 2 | no | no |
| Br1898 | caud-dlpf | either | 1 of 2 | no | no |
| **Br1249** | dlpf-hipp | either | 1 of 2 | **no** | **yes** |
| **Br1693** | caud-dlpf | either | 1 of 2 | **no** | **yes** |
| **Br1700** | dlpf-hipp | either | 1 of 2 | **no** | **yes** |

**This inverts the concern T4 was written about.**

- **No single swapped tube explains any of the three readmitted donors.** Each is
  flagged in exactly one pair while the other pair sharing each candidate region
  is clean — and for Br1249 and Br1700 that clean pair is **rank 1 of ~102**, the
  strongest possible self-match. A DLPFC swap for Br1249 would have to break
  caudate-DLPFC (rank 4, z +1.567); a hippocampus swap would have to break
  caudate-hippocampus (rank 1, z +4.820). Neither happens. The flag is
  incompatible with the mechanism it is screening for.
- **The only two donors whose pattern IS swap-coherent are Br1134 (hippocampus)
  and Br1459 (DLPFC), and neither is a readmitted donor.** If anything in this
  screen deserves follow-up, it is those two — and they have nothing to do with
  the blacklist.

Supporting detail: the flags are not uniformly distributed across pairs. DLPFC-
hippocampus carries 6 of 11 flags at a 5.2% rate, against 2.9% for caudate-DLPFC
and 2.0% for caudate-hippocampus, and it also has the lowest median `z_self`
(+2.57 vs +3.12 and +3.49). DLPFC and hippocampus are the two regions sharing 115
of 118 donors and sitting in the same sequencing batches, so the pair with the
most candidates and the least region contrast is also the noisiest. Five of the
nine flagged donors are flagged **only** there.

Scale: `z_self < 0` occurs in 11 of 320 rows (3.4%). The 5th percentile of
`z_self` is +0.183 and the median is +2.915, so the flag threshold sits inside the
lower tail of a broad distribution rather than at a natural break.

## 5. What would actually settle it

The screen names its own remedy: **the unmasked C->T SNP CpGs**. Those are the
genotype-driven sites Module 01 removes, and they are where a methylation
fingerprint has real power. A check on them would be a verdict rather than a
nomination. It needs the pre-mask CpG matrix, so it is a new read of the source
BSobjs rather than anything derivable from a sealed run — a bounded piece of work,
and the only thing that converts this from a judgement call into a measurement.

## 6. Recommendation

**Clear all eight and affirm `config/cohorts.yml:300`, on the record, with the
reasoning in section 4.** Grounds:

1. Only three of the eight are flagged at all, and `Br1927` — named in T4 — is
   clean in all three pairs and is its own best match out of 103.
2. None of the three flags is coherent with a sample swap. Each breaks one pair
   while its companion pair, sharing the same candidate tube, is clean and in two
   cases is rank 1.
3. The two swap-coherent donors in the cohort are not readmitted donors.
4. The association rests on three events; the donor x region-pair interval
   includes 1, and neither unit was prespecified.
5. The screen's author states it cannot adjudicate an individual donor.
6. The original grounds for readmission are unrelated to QC: the blacklists
   existed to reconcile a stale AA phenotype file, v2 reads one phenotype table,
   and all eight donors have usable WGBS and genotypes.

**Two things to record alongside the decision, because they are real findings and
would otherwise be lost:**

- **Br1134 and Br1459 are the screen's actual nominations** — coherent with a
  hippocampus and a DLPFC swap respectively. They are not blacklist donors and
  have never been assessed. Either file them for the C->T follow-up or record
  that the screen is too weak to act on them; do not leave them unmentioned.
- **T4's numbers should be retired, not corrected in place.** They describe a
  deleted run. Replace them with section 3's two tables and the unit caveat.

If instead any donor is to be excluded, note that it re-derives Module 01 and
therefore `vmr_set_id`, and every module reruns. On this evidence that cost is not
purchased by a corresponding gain in confidence.

---

## Decision as taken

**2026-09-25, Kynon J.M. Benjamin (PI): all eight donors cleared;
`config/cohorts.yml` affirmed.**

`Br1249 Br1303 Br1371 Br1552 Br1693 Br1700 Br1883 Br1927` remain INCLUDED in v2.
`sample_blacklist` stays null for all three regions. The assertion that the legacy
blacklists "were never a QC exclusion" stands, now with the cross-region integrity
screen read against it rather than left unaddressed.

### What this closes

- **T4 is closed.** Its recorded p = 0.0064 / OR 9.8 describe the deleted run
  `fig-all-20260826`/`-a` and are retired, not corrected; section 3 replaces them.
- **No rerun follows.** The donor set is unchanged, so Module 01 does not
  re-derive, `vmr_set_id` is stable in all three regions, and no downstream module
  is invalidated. This removes the last item that could have invalidated the
  Module 05 runs awaiting acceptance.

### What this does NOT close, and must not be lost

- **`Br1134` and `Br1459` are the screen's only swap-coherent nominations** —
  consistent with a hippocampus and a DLPFC mislabel respectively (section 4).
  Neither is a readmitted donor and neither has been assessed. They are **open**,
  not cleared by this decision, which was scoped to the eight blacklist donors.
- **The C->T follow-up remains the only test that could adjudicate an individual
  donor** (section 5). Until it is run, this screen nominates and does not convict
  — for `Br1134` and `Br1459` as much as for the eight.

Neither of those two points is a blocker: the screen is a supplementary QC figure
and nothing downstream consumes it (`_h/05_qc_sample_integrity.R:30`).
