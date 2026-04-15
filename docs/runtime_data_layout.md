# Runtime data layout

The package resolves data from package-managed roots under `data/` and from the configured `references/` and `public/` trees.

## Canonical runtime roots

- `data/tables/`
  - shipped CSV defaults such as `table_coeff.csv`, `table_filter.csv`, `table_meta.csv`, `table_param.csv`
- `data/catalog/`
  - operator-maintained registries such as `list_deployment.csv` and `list_deployment_hr.csv`
- `data/config_files/`
  - deployment and platform JSON metadata, including patch files
- `data/raw_smru_data_odv/`
  - low-resolution ODV zip/text inputs
- `data/raw_smru_hr_data/`
  - high-resolution raw CTD inputs used by the FR branch
- `references/`
  - reference datasets such as WOD, CORA, regression fixtures, and the last stable MEOP release
- `final_dataset_prof/`
  - generated `lr0`, `lr1`, `hr0`, `hr1`, `hr2`, `fr0`, `fr1`, and `traj` products
- `plots/`
  - diagnostics figures
- `data/batch/`
  - resumable batch state, logs, and run summaries

## Low-resolution raw data

The package accepts any of these staged directly under `data/raw_smru_data_odv/`:

- `<deployment>_ODV.zip`
- `<deployment>_ODV.txt`
- `<deployment>_CTD_ODV.txt`
- `<deployment>_FL_ODV.txt`

If a zip archive is present, the import step extracts it in place.

## Expected high-resolution raw data

High-resolution CTD files are resolved from `list_deployment_hr.csv` and expected at:

- `data/raw_smru_hr_data/<year>/<optional_prefix><instr_id>_ctd.txt`

The helper `meop_process.io.hr_ctd.resolve_hr_ctd_path()` resolves that location.

## Current processing ownership

Python-owned today:

- deployment and tag discovery
- ODV import and profile indexing
- `lr0` generation
- QC/filtering
- location-adjustment placeholder
- `hr0`, `hr1`, `lr1`
- `apply_adjustments`
- `fr0`, `fr1`
- `hr2`
- standard diagnostics figures
- resumable batch processing and summary-table refresh

Still deferred or intended for later redesign:

- specialized analyst comparison figures
- richer diagnostic and analyst-comparison workflows
- publication/reporting helpers that need a modern Python redesign
