# MEOP data layout

The package resolves data through explicit runtime roots owned by the Python package.
No root-level mirrors are required.

## Data classes and canonical locations

1. **Packaged CSV defaults**
   - shipped inside `src/meop_process/resources/tables/`
   - synchronized at runtime into `data/tables/`

2. **Operator-maintained catalogs**
   - `list_deployment.csv`
   - `list_deployment_hr.csv`
   - optional `list_wmo_numbers.csv`
   - canonical runtime location: `data/catalog/`

3. **Runtime / private process data**
   - low-resolution raw ODV archives and text files under `data/raw_smru_data_odv/`
   - high-resolution CTD text files under `data/raw_smru_hr_data/`
   - deployment/platform JSON metadata under `data/config_files/`
   - crawl and CLS location updates under `data/crawl_locations/` and `data/smooth_cls_locations/`

4. **Reference datasets**
   - large or private datasets under `references/`
   - expected fixed subtrees include WOD, CORA, and the last stable MEOP release

5. **Outputs**
   - processed netCDF under `final_dataset_prof/`
   - diagnostics under `plots/`
   - maps under `maps/`
   - LaTeX products under `doc_latex/`
   - public release tree under `public/<version>/`

## Runtime resolution order

### Tables

1. `data/tables/<name>.csv`
2. packaged defaults from `src/meop_process/resources/tables/`

### Catalogs

1. `data/catalog/<name>.csv`
2. when loading deployments specifically, missing or incomplete information can be filled from
   JSON files in `data/config_files/`

### Raw low-resolution data

1. `data/raw_smru_data_odv/<deployment>_ODV.zip`
2. `data/raw_smru_data_odv/<deployment>_ODV.txt`
3. `data/raw_smru_data_odv/<deployment>_CTD_ODV.txt`
4. `data/raw_smru_data_odv/<deployment>_FL_ODV.txt`

### Raw high-resolution data

1. `data/raw_smru_hr_data/<year>/<prefix><instr_id>_ctd.txt`

## Current workflow boundary

Python now owns the complete operational workflow currently implemented in the package:

- deployment and HR catalog loading
- JSON metadata loading from `data/config_files/`
- low-resolution raw ODV import and parsing
- `lr0`, `lr1`, `hr0`, `hr1`, `hr2`, `fr0`, `fr1`, and `traj`
- delayed-mode adjustments
- standard diagnostics
- resumable batch reruns
- incremental refresh of `list_tags.csv` and `list_deployments.csv`

Still deferred for later redesign:

- specialized analyst comparison figures
- richer analyst comparison tools
- publication/reporting helpers such as ODV export and legacy LaTeX report generation
