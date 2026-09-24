# Decision pilot: the executed Module 05 covariate model vs the locked M3a

**Status: evidence for a PI decision. Not a fix, not a production run.**
Nothing under `_m/` was created, modified or re-sealed. No run ID was minted.
All outputs live in scratch; the harness is `_h/pilot_covariate_model.{py,sh}`
and `tests/pilot_signal_check.py`.

## The question

`config/covariates.yml:primary_meqtl` is locked (`lock_status: locked_M3a`,
`lock_date: 2026-08-01`) to

```text
M3a = agedeath + sex + primarydx + snpPC1-5 + methPC1-5
```

`lock_rationale`: *"M3a improves λ_NS vs M0 and increases FDR-significant CpG
discoveries in all three regions while retaining external meQTL enrichment."*

`_h/01b_prepare_meqtl_inputs.py:145` sets `n_pc = 3` and never adds a
methylation PC, so what the three accepted runs actually fitted is

```text
executed = agedeath + sex + primarydx + snpPC1-3
```

Nothing reads `config/covariates.yml` in Module 05 or in `00_shared/`. The only
reference in the repository is `00_shared/runid.R:73`, which writes
`config_covariates_sha256` into the manifest. The lock is therefore attested by
every sealed run and enforced by nothing, which is how a covariate model drifted
from a PI decision without a gate firing. AGENTS.md §7.5 makes genomic inflation
a gate and §12 makes the covariate model a PI decision, so this pilot supplies
the number and does not make the choice.

## The chromosome, and how far it generalises

**chr10.** Prespecified grounds: a mid-sized autosome with the fourth-largest
tested-CpG count in the caudate run, no MHC, far from the small-chromosome
regime where the distal-null tail thins. Chosen before any arm was fitted, and
not by its λ.

| quantity | chr10 | genome-wide (caudate) | share |
|---|---:|---:|---:|
| member CpGs | 10,164 | 191,946 | 5.30% |
| VMRs with a tested CpG | 604 | 11,352 | 5.32% |
| nominal cis pairs | 30,034,746 | — | — |
| distal (> 400 kb) pairs | 5,758,046 | — | — |

λ is not uniform across chromosomes, which bounds the extrapolation of the
*level*. Recomputed from the sealed caudate run's own nominal parquets:

| chromosome | λ distal (> 400 kb) | λ all pairs |
|---|---:|---:|
| chr10 | 1.1651 | 1.3941 |
| chr11 | 1.1081 | 1.3302 |
| chr12 | 1.1713 | 1.4175 |
| chr22 | 1.0873 | 1.2337 |
| **genome-wide pooled (the reported gate value)** | **1.1391** | **1.3762** |

A single chromosome therefore pins the *difference between arms* far better than
it pins the level. chr10's executed level is above the pooled 1.1391 and above
`genomic_inflation.max: 1.15`; that is between-chromosome variation, and it is
why the gate is applied to the pooled figure rather than per chromosome.

## The comparison

All five arms share the sealed run's tested-CpG set and phenotype BED, so only
the covariate design and the genotype-QC donor set move. λ is the module's own
statistic — `_h/04_qc_plots.py::genomic_inflation` over nominal pairs with
|CpG-to-SNP distance| > `genomic_inflation.distal_min_distance_bp` (400 kb) —
imported from that file rather than reimplemented. It is computed over **all**
of chr10's pairs rather than the module's 5 M subsample, which is the same
estimand with less noise.

| arm | covariates | n | λ distal | λ all pairs | tested CpGs | nominal pairs | distal pairs | Storey q≤0.05 (chr10 family) | BH q≤0.05 | π0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **executed** (snpPC1-3) | 6 | 153 | **1.1651** | 1.3941 | 10,164 | 30,034,746 | 5,758,046 | 7,197 | 4,533 | 0.137 |
| m0 (snpPC1-5) | 8 | 153 | 1.1701 | 1.3931 | 10,164 | 30,034,746 | 5,758,046 | 6,111 | 4,482 | 0.216 |
| **m3a** (snpPC1-5 + methPC1-5) | 13 | 153 | **1.1412** | 1.3749 | 10,164 | 30,034,746 | 5,758,046 | 6,260 | 4,921 | 0.226 |
| executed_keep | 6 | 153 | 1.1619 | 1.3889 | 10,164 | 30,969,219 | 5,947,863 | 6,625 | 4,516 | 0.167 |
| m3a_keep | 13 | 153 | 1.1425 | 1.3697 | 10,164 | 30,969,219 | 5,947,863 | 6,218 | 4,903 | 0.225 |

Donor n is **153 in every arm**. The donor set was never the thing at risk: the
mapping stage intersects phenotype columns with pgen samples, so the regression
always used the 153 caudate donors. What the missing `--keep` changed is which
*SNPs* exist to be tested.

Reading the covariate axis:

- **M3a lowers λ by 0.0239** (1.1651 → 1.1412, −2.1%). Direction agrees with the
  lock's rationale. The 2026-08-01 lock recorded Δλ_NS of −0.009 (caudate),
  −0.035 (DLPFC), −0.041 (hippocampus) on a *different* statistic (λ over
  non-significant lead p-values), so −0.024 on the distal-null statistic is the
  same order of magnitude.
- **The improvement is entirely the latent factors, not the extra ancestry PCs.**
  m0 (snpPC1-5, no methPC) gives λ = 1.1701, marginally *worse* than the executed
  snpPC1-3. snpPC4-5 contribute nothing; methPC1-5 move λ from 1.1701 to 1.1412.
  This is why the m0 arm was added, and it is the reason a two-arm pilot would
  have mis-attributed the effect.

## Did methylation PCs buy λ by eating cis signal?

No, on this chromosome. The one number that says otherwise is the least
trustworthy in the table.

The Storey count **fell** under M3a (7,197 → 6,260, −13.0%) while the BH count
**rose** (4,533 → 4,921, +8.6%). The two disagree in direction because π0 is
estimated from the p-value tail and moves from 0.137 to 0.226 between the arms.
`fdr_family: per_brain_region` means the module's real family spans all 22
autosomes; a one-chromosome Storey family on 10,164 CpGs is **not** the module's
family, and π0 is the least stable quantity in it. So the Storey column here
should not be read as a discovery count.

Everything that does not depend on π0 says M3a gains:

```
CpGs below a FIXED permutation-p threshold (no pi0, no FDR)
arm            p<=1e-2  p<=1e-3  p<=1e-4  p<=1e-5  p<=1e-8
executed          4230     3756     3311     2980     2567
m0                4238     3737     3319     2975     2538
m3a               4519     3865     3481     3148     2591
executed_keep     4249     3732     3301     2969     2559
m3a_keep          4503     3885     3462     3131     2583
```

M3a is ahead at every threshold (+6.8%, +2.9%, +5.1%, +5.6%, +0.9%), and the BH
set movement is +557 gained against 169 lost. Rank agreement between the two
models' permutation p-values is Spearman ρ = 0.925 over all 10,164 CpGs.

There *is* a real degrees-of-freedom cost, and it is worth stating because it is
the mechanism a reviewer would ask about. Among the CpGs the executed model
already called at BH 0.05, evidence is slightly weaker under M3a (median change
−0.146 in −log10 p; 54.3% of them weaker). But m0 is worse on the same measure
(−0.226; 77.2% weaker) while gaining nothing, so most of that attenuation is the
price of extra covariates in 153 donors rather than anything specific to the
latent factors. Net of both effects M3a still discovers more.

## What the methylation PCs actually are

`config/covariates.yml` names **no file** for methPC1-5. It names a method:
`latent_factor_policy: "Locked primary uses methPC1–5 (PCA on M0-residualized
CpG phenotypes)"`. That is a genuine underspecification and a finding in its own
right — the lock does not fix the CpG set, the subsample size, the seed, the
standardization, or whether the M0 residualization uses snpPC1-5.

The pilot resolved it by reimplementing the v1 implementation of that method,
`meqtl-validation/01_cpg_meqtl_mapping/_h/09_estimate_latent_factors.py`, under
its own seed (20260730) and 50,000-CpG subsample, applied to the **v2** tested-CpG
phenotypes pooled over all 22 autosomes (191,946 CpGs × 153 donors). The v1
factor tables were deliberately not reused: they were estimated on the v1 VMR
catalog's CpGs. methPC1-5 explain 15.8% of residual CpG variance
(5.96/4.36/2.37/1.76/1.32% individually).

**A consideration the finding as filed does not mention, and which may matter
more than λ.** Regressing each factor on the region's RNA MuSiC proportions:

| factor | R² on RNA MuSiC proportions | strongest single correlation |
|---|---:|---|
| methPC1 | **0.721** | Oligo ρ = +0.766 (p = 9.2e-31); D1-SPN ρ = −0.63; D2-SPN ρ = −0.62 |
| methPC2 | 0.062 | OPC ρ = +0.218 |
| methPC3 | 0.191 | Astro ρ = −0.325 |
| methPC4 | 0.071 | Inhib ρ = −0.150 |
| methPC5 | 0.076 | OPC ρ = −0.185 |

methPC1 is largely collinear with estimated cell composition. This is a
bulk-tissue correlation between a methylation PC and an RNA-derived proportion
estimate; it establishes collinearity and **not** a cell type of origin, which
§2.3 forbids inferring. But it means adopting M3a imports a substantial
cell-composition adjustment into the **primary** meQTL model, and
`config/covariates.yml` currently says `cell_composition: sensitivity_only`,
with the cell-adjusted designs registered as the separate, gated M5 and M6d
sensitivities. It would also weaken the M6d contrast, since M6d is defined as
`M3a + dnamCellPC1-3` and much of what those PCs adjust for would already be in
M3a. That interaction is a PI matter, not a statistical one.

## The discarded `--keep`, and a second defect inside it

`_h/01b_prepare_meqtl_inputs.py:203-204` writes `{chrom}.keep` from the analysis
donors; the plink2 call at 206-214 passes `--pfile --chr --maf --geno --hwe
--threads --make-pgen --out` and never `--keep`. The sealed chr10 log confirms
the consequence directly:

```
526 samples (...) loaded from .../chr10_src.psam
--geno: 2596 variants removed due to missing genotype data.
346532 variants removed due to allele frequency threshold(s)
416659 variants remaining after main filters.
```

The intended set is the 153 caudate donors (`manifest.tsv:n_donors = 153`). The
source `TOPMed_LIBD.AA` pfile holds **526** AA donors, so MAF, missingness and
HWE were evaluated over **373 donors outside the estimation set**.

`00_shared/locus_io.R:33-34` documents the opposite convention for Modules 02
and 03 — the filters are computed *after* the group restriction, deliberately, so
that they are computed within the estimation group. Module 05 is the divergence.

**The keep file is also malformed, and that is why the run completed.** It writes
`{donor}\t{donor}` — the FID twice — but the psam's IID is a chip barcode
(`Br2585  3998646007_R01C01`). plink2 reads a two-column `--keep` as FID/IID, so
it matches nothing:

```
$ plink2 --pfile .../chr10_src --keep .../chr10.keep --chr 10 --write-samples
--keep: 0 samples remaining.
Error: No samples remaining after main filters.
```

Had `--keep` been passed as written, every array task would have failed. A usable
restriction has to pair each FID with its psam IID; the pilot does that, and
plink2 then reports `--keep: 153 samples remaining`.

Effect on chr10's variant set, correct (153-donor) filters versus executed
(526-donor) filters:

| | variants |
|---|---:|
| executed pfile (no `--keep`) | 416,659 |
| correct pfile (`--keep`, 153 donors) | 414,363 |
| shared | 398,827 |
| **over-included** — selected using donors outside the estimation set | **17,832** |
| **wrongly excluded** — pass the locked QC in the 153, never tested | **15,536** |

Decomposing each side:

- Of the 17,832 over-included, 16,085 (90.2%) have MAF < 0.05 in the 153 and are
  removed again by tensorqtl's in-sample MAF filter. The remaining **1,747 are
  actually tested, and all 1,747 exceed the locked `missingness_max: 0.05` in the
  153 donors** — tensorqtl re-applies MAF in sample but never missingness, so
  nothing downstream catches them.
- Of the 15,536 wrongly excluded, 14,958 failed `--maf` in the 526 and 578 failed
  `--geno` there. All of them pass the locked QC in the 153.

Net: the executed run tested about **400,574** SNPs on chr10 where the locked QC
gives **414,363** — roughly **3.4% fewer testable SNPs**, plus 1,747 tested SNPs
that violate the configured missingness threshold. The pilot's pair counts agree:
30,034,746 nominal pairs without `--keep` against 30,969,219 with it (+3.1%),
*more* pairs from *fewer* pfile variants, exactly as the in-sample MAF filter
predicts.

**But it does not move λ.** executed 1.1651 → executed_keep 1.1619 (−0.0032);
m3a 1.1412 → m3a_keep 1.1425 (+0.0013). Both are far smaller than the covariate
effect and far smaller than between-chromosome variation. BH counts move by
≤ 0.4% (4,533 → 4,516; 4,921 → 4,903). So this defect is a **denominator and
QC-compliance** problem — which SNPs were eligible, and 1,747 that should not
have been — rather than a calibration problem. It is not a reason to rerun on its
own, and it is a reason that any rerun undertaken for another purpose should
carry the fix.

## The three-way reading

The pilot measures one chromosome. Stating all three branches, as asked, and
without recommending one:

**If M3a is materially better.** It is better here, but by 0.024 in λ on 5.3% of
the CpGs — "material" is the PI's threshold, not the pilot's. Taken as material,
the lock stands, Module 05 reruns under M3a, and λ must be re-reported (chr10's
−0.024 would put a pooled λ near 1.115 if it carried, i.e. further inside the
0.90–1.15 band, but that is an extrapolation from one chromosome and the other
two regions are untested). The rerun propagates: `07/_h/03_test_coupling.R:74`
derives `any_meqtl_support` from Module 05's per-VMR counts, so 07 reruns, then
08 (`01_cross_region_replication.R:146`), then 09 (`03_risk_variant_cpg_meqtl.py`
and `08_prepare_coloc_regions.py` both read the burden run), then 11
(`08_figure4_meqtl_coupling.R`). Five stages. Per TASKS.md this is the single
item that can lengthen the critical path, and the F3 PSI re-keying should ride
the same Module 07 rerun rather than triggering a second one.

**If the executed model is materially better.** Nothing in these numbers supports
that: λ is worse (1.1651 vs 1.1412) and discoveries are fewer at every π0-free
threshold. If the PI nonetheless ratifies the executed model — on parsimony, on
keeping cell composition out of the primary per `cell_composition:
sensitivity_only`, or on rerun cost — then **`config/covariates.yml` has to be
amended to say so**, because the executed design is not a named model in that
file. `sensitivity_models` has no snpPC1-3 entry, and `primary_meqtl.ancestry_pcs`
lists snpPC1-5, so "the executed model" is currently an unnamed design that
exists only as the literal `3` at `01b:145`. No Module 05 rerun for the covariate
question; 07/08/09/11 are untouched; the Methods text changes and the lock's
rationale is superseded. A gate that compares the executed design against the
config should be added either way, so this cannot recur silently.

**If there is no material difference.** Then λ does not decide it, and the
deciding consideration is **which cell-composition strategy the primary meQTL
model is supposed to have**. methPC1 is 72% explained by RNA MuSiC proportions,
so the choice between these two models is substantively a choice about whether
the primary scan is cell-adjusted — which `config/covariates.yml` presently
assigns to the gated M5/M6d sensitivities and AGENTS.md §7.4/§7.5 treat as a
locked sensitivity rather than part of the primary. A secondary tiebreak is the
record: §9 requires a run to carry its configuration, and a manifest that
checksums a lock the run did not follow does not satisfy that, whichever way the
PI resolves it.

## What this pilot does not resolve

- One chromosome, 5.3% of CpGs, caudate only. DLPFC and hippocampus are untested,
  and the lock's original evidence was strongest in those two regions
  (Δλ_NS −0.035 and −0.041 against caudate's −0.009).
- The module's FDR family cannot be reproduced on one chromosome, so the
  significant-CpG *count* under the module's own rule is out of reach here. The
  π0-free thresholds and BH are the substitutes, and they agree with each other.
- The effect on the module's actual endpoint — `proportion_cpgs_with_sig_meqtl`
  per VMR, and the `local_snp_contribution_score_z` coefficient — is not
  estimated. chr10 carries 604 VMRs and the burden model is fitted genome-wide.
- External meQTL enrichment, the lock's third criterion, was not recomputed.

## Reproducing

```bash
R=/gpfs/.../dna-methylation-heritability          # or a worktree
S=/scratch/.../pilot
conda run -p /projects/p32505/opt/envs/genomics python \
    $R/05_cpg_meqtl_burden/_h/pilot_covariate_model.py methpcs --out $S
sbatch --chdir=$S/slurm \
    --export=ALL,REPO=$R,PILOT_OUT=$S,CHROM=10 \
    $R/05_cpg_meqtl_burden/_h/pilot_covariate_model.sh
conda run -p /projects/p32505/opt/envs/genomics python \
    $R/05_cpg_meqtl_burden/tests/pilot_signal_check.py
```

Jobs actually submitted: `7172434` (arms executed, m3a, executed_keep, m3a_keep;
1 h 47 m 45 s) and `7172553` (arm m0; 30 m 12 s), both `genomics` partition,
8 CPUs, 64 GB, `COMPLETED`. The `executed` arm is a reproduction check and it
passes twice over: its covariate file is byte-identical to the sealed run's
(md5 `61adc2091d03bf93fefd34ef11dd76d9`), and its λ matches the value recomputed
independently from the sealed nominal parquet to four decimal places
(1.1651 both ways).
