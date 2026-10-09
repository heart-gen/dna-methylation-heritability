# zenodo — the per-module deposit

`zenodo_manifest.tsv` lists one zip per module, built by
`supplementary_data/_h/build_release.py --zip`. The zips themselves are
gitignored; upload them as a new version of the Zenodo record
[10.5281/zenodo.20547606](https://doi.org/10.5281/zenodo.20547606).

Each zip holds, with paths relative to the repository root:

- every **manuscript file** of the module (also tracked in this repository, in
  git or Git LFS), so the DOI is self-contained; and
- the module's **reproducibility extras** that no manuscript text draws on
  (run provenance, uncited per-VMR tables, annotations, LD scores). These are
  on Zenodo only.

No zip contains individual-level methylation, genotype or covariate data
(dbGaP phs000979), SVG copies of figures, logs, or bulk intermediates. The
file-by-file assignment, with the reason for each file's tier, is
`supplementary_data/release_manifest.tsv`.

| column | meaning |
|---|---|
| `archive` | zip file name |
| `module` | the module it packages |
| `bytes`, `sha256` | of the zip, for verifying the upload |
| `n_members` | files in the zip |
| `n_manuscript_members` | of those, files also tracked in git/LFS |
| `sd_numbers` | Supplementary Data items the zip contains |
| `runs` | accepted run IDs the files come from |

Rebuild after any acceptance changes:

```bash
python3 supplementary_data/_h/build_release.py --zip --stage --check
```
