"""Read configuration from a run's immutable snapshot.

The submit drivers snapshot _h/ into runs/{RUN_ID}/code/_h so an edit to _h/
cannot change a queued run. Configuration was not covered by that: stages
re-read config/ from the live working tree, so a config edit -- or a git branch
switch that removed the file, as happened on 2026-09-02 -- changed or broke a
run already in flight. The driver now snapshots config/ alongside the code, and
this loader prefers the snapshot, falling back to the live file for a hand-run
stage.
"""
from __future__ import annotations

from pathlib import Path

import yaml


def load_run_config(name: str, run_dir: Path, root: Path) -> dict:
    snap = Path(run_dir) / "code" / "config" / f"{name}.yml"
    src = snap if snap.exists() else Path(root) / "config" / f"{name}.yml"
    if not src.exists():
        raise SystemExit(f"Config not found: {src}")
    return yaml.safe_load(src.read_text())
