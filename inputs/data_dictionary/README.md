# Data dictionary

Phase 0 formal data audit for strengthening our findings.

## Run

Run from `../_m/`:

```bash
mkdir -p logs
sbatch ../_h/step_1.sh
```

Or locally:

```bash
python3 ../_h/00_build_data_dictionary.py
```

Reads `config/paths.yml` when PyYAML is available.

## Outputs (`_m/`)

| File | Description |
|---|---|
| `sample_overlap.tsv` | Phenotype × genotype × cell × CpG matrix overlap by race/region |
| `covariate_dictionary.tsv` | Covariate roles and missingness policy |
| `genome_build_audit.tsv` | Genome build and liftOver needs per dataset |
| `public_meqtl_resources.tsv` | External brain meQTL resource catalog |
| `region_input_summary.tsv` | VMR / predictability / CpG inventory |
| `data_dictionary_reproducibility.tsv` | Execution metadata |

See also `_h/README.md` for download notes.
