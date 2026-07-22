# MEOP data layout

The package resolves data through explicit runtime roots owned by the Python package.
No root-level mirrors are required.

## Data classes and canonical locations

1. **Packaged CSV defaults**
   - shipped inside `src/meop_process/_bundled/tables/`
   - synchronized at runtime into `data/tables/`
   - `table_coeff_comment_codes.csv` documents the standardized `table_coeff.csv` comment codes

2. **Operator-maintained catalogs**
   - `list_deployment.csv`
   - `list_deployment_hr.csv`
   - optional `list_wmo_numbers.csv`
   - canonical runtime location: `data/catalog/`

3. **Runtime / private process data**
   - low-resolution raw ODV archives and text files under `data/data_raw/raw_smru_data_odv/`
   - high-resolution CTD text files under `data/data_raw/raw_smru_hr_data/`
   - deployment/platform JSON metadata under `data/data_raw/config_files/`
   - crawl and CLS location updates under `data/data_raw/crawl_locations/` and `data/data_raw/smooth_cls_locations/`

4. **Outputs**
   - processed profile netCDF under `data/data_prof/`
   - processed trajectory netCDF under `data/data_traj/`
   - diagnostics under `data/plots_by_tags/`
   - deployment and overview plots under `data/plots_by_deployments/` and `data/plots_overview/`
   - maps under `data/maps/`
   - public release tree under `public/<version>/`

There is no longer any transition fallback to removed root-level output directories such as `final_dataset_prof/` or `plots/`.

## Runtime resolution order

### Tables

1. `data/tables/<name>.csv`
2. packaged defaults from `src/meop_process/resources/tables/`

### Catalogs

1. `data/catalog/<name>.csv`
2. when loading deployments specifically, missing or incomplete information can be filled from
   JSON files in `data/data_raw/config_files/`

### Raw low-resolution data

1. `data/data_raw/raw_smru_data_odv/<deployment>_ODV.zip`
2. `data/data_raw/raw_smru_data_odv/<deployment>_ODV.txt`
3. `data/data_raw/raw_smru_data_odv/<deployment>_CTD_ODV.txt`
4. `data/data_raw/raw_smru_data_odv/<deployment>_FL_ODV.txt`

### Raw high-resolution data

1. `data/data_raw/raw_smru_hr_data/<year>/<prefix><instr_id>_ctd.txt`

## Current workflow boundary

Python now owns the complete operational workflow currently implemented in the package:

- deployment and HR catalog loading
- JSON metadata loading from `data/data_raw/config_files/`
- low-resolution raw ODV import and parsing
- geographic location adjustment
- `lr0`, `lr1`, `hr0`, `hr1`, `hr2`, `fr0`, `fr1`, and `traj`
- delayed-mode adjustments
- standard diagnostics
- resumable batch reruns
- incremental refresh of `list_tags.csv` and `list_deployments.csv`

By default, successful production processing retains final profile products (`hr1`, `lr1`, and `hr2`) and prunes rebuildable intermediates (`lr0`, `hr0`, `fr0`, and `fr1`).
Set `defaults.processing.keep_intermediate` or pass `--keep-intermediate-products` for calibration/debug runs.

Batch reruns use `data/batch/latest/deployment_status.json` as persistent state, but successful entries are pruned automatically if their canonical outputs under `data/data_prof/` have disappeared.

Still deferred for later redesign:

- specialized analyst comparison figures
- richer analyst comparison tools
- publication/reporting helpers such as ODV export and legacy LaTeX report generation
