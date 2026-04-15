# MEOP package structure

The repository now uses a `src/`-layout Python package with a single pure-Python execution path.

## Current package shape

Main modules:

- `src/meop_process/api.py`
- `src/meop_process/cli.py`
- `src/meop_process/catalog/`
- `src/meop_process/io/`
- `src/meop_process/processing/`
- `src/meop_process/plotting/`
- `src/meop_process/batch/`
- `src/meop_process/metadata/`
- `src/meop_process/data/`

## Main public workflow functions

- `process_tags`
- `create_fr0`
- `create_hr2`
- `apply_adjustments`
- `generate_diagnostics`
- `run_all_deployments`
- `update_metadata_summaries`

## Design choices in the cleaned package

- a single Python execution path;
- no dependency on external language runtimes;
- no root-level table or catalog mirrors;
- all managed runtime data live under `data/`;
- legacy specialized plotting entrypoints removed from the main package workflow.
