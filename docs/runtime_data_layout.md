# Runtime data layout

The package resolves runtime data from package-managed roots under `data/` and from the configured `public/` tree.

## Canonical runtime roots

- `data/tables/`
  - shipped CSV defaults such as `table_coeff.csv`, `table_filter.csv`, `table_meta.csv`, `table_param.csv`
- `data/catalog/`
  - operator-maintained registries such as `list_deployment.csv` and `list_deployment_hr.csv`
- `data/data_raw/config_files/`
  - deployment and platform JSON metadata, including patch files
- `data/data_raw/raw_smru_data_odv/`
  - low-resolution ODV zip/text inputs
- `data/data_raw/raw_smru_hr_data/`
  - high-resolution raw CTD inputs used by the FR branch
- `data/data_prof/`
  - generated `lr0`, `lr1`, `hr0`, `hr1`, `hr2`, `fr0`, and `fr1` profile products
- `data/data_traj/`
  - generated trajectory products
- `data/plots_by_tags/`
  - per-tag diagnostics figures grouped by deployment
- `data/plots_by_deployments/`
  - deployment-level recap plots and summaries
- `data/plots_overview/`
  - cross-deployment overview summary plots
- `data/maps/`
  - map outputs
- `data/batch/`
  - resumable batch state, logs, and run summaries

## Optional external reference data

- `cora_dir` (set in `configs.json` under `defaults.references.cora_dir` or `configs.<machine>.references.cora_dir`)
  - directory of CORA 10°×10° tiled NetCDF files, e.g.
    `/path/to/CORA_ncfiles/CORA_lon040W_lat70S.nc`
  - used by `meop-compare --plot1` to generate CORA-based T/S calibration plots
  - not required for the core processing pipeline

## Low-resolution raw data

The package accepts any of these staged directly under `data/data_raw/raw_smru_data_odv/`:

- `<deployment>_ODV.zip`
- `<deployment>_ODV.txt`
- `<deployment>_CTD_ODV.txt`
- `<deployment>_FL_ODV.txt`

If a zip archive is present, the import step extracts it in place.

## Expected high-resolution raw data

High-resolution CTD files are resolved from `list_deployment_hr.csv` and expected at:

- `data/data_raw/raw_smru_hr_data/<year>/<optional_prefix><instr_id>_ctd.txt`

The helper `meop_process.io.hr_ctd.resolve_hr_ctd_path()` resolves that location.

## Current processing ownership

Python-owned today:

- deployment and tag discovery
- ODV import and profile indexing
- `lr0` generation
- QC/filtering
- geographic location adjustment from crawl, CLS, and SMRU sources
- `hr0`, `hr1`, `lr1`
- `apply_adjustments`
- `fr0`, `fr1`
- `hr2`
- standard diagnostics figures
- resumable batch processing and summary-table refresh

Batch state is stored in `data/batch/latest/deployment_status.json` and reconciled against the canonical output tree at batch startup. Successful entries whose outputs have been deleted are dropped from state before skip decisions are made.

Still deferred or intended for later redesign:

- publication/reporting helpers that need a modern Python redesign
