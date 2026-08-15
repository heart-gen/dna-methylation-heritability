# Manuscript Supplementary Data builder

`_h/build_supplementary_data.py` creates Supplementary Data 10–15 from the
final outputs of the meQTL validation workflow.

Run from the repository root:

```bash
python3 meqtl-validation/12_supplementary_data/_h/build_supplementary_data.py
```

For the SLURM repair rerun, submit `_h/step_1.sh` only after every schema-v2
analysis dependency succeeds.

By default, generated archives and their release manifest are written beneath
`meqtl-validation/12_supplementary_data/_m/`, which is intentionally ignored by
Git. Use `--output PATH` to stage them elsewhere.

The builder uses a fixed archive timestamp and records file sizes, SHA-256
checksums, schemas, and row counts. It rejects duplicate table headers,
individual-level identifier columns, and cluster-local filesystem paths. It
does not package logs, archived model runs, sample/donor lists, internal
strategy or readiness decisions, individual-level molecular matrices, or raw
third-party meQTL/GWAS resources.
