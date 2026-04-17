# Runtime resolution and task ownership

This document describes how the cleaned package resolves data and which steps are currently implemented.

## `process_tags()` now runs in five Python stages

1. **Runtime preparation**
   - ensure runtime directories exist
   - optionally synchronize deployment/platform JSON files into `data/data_raw/config_files/`

2. **Selection and raw-data discovery**
   - resolve deployment metadata with `load_info_deployment()`
   - import or extract low-resolution ODV files from `data/data_raw/raw_smru_data_odv/`

3. **Cleanup and metadata preparation**
   - remove prior outputs for the target deployment or tag
   - patch processed-file metadata from `table_meta.csv`

4. **Core processing**
   - `create_ncargo` -> `lr0`
   - QC/filtering
   - geographic location adjustment
   - `create_hr0`
   - `create_fr0`
   - `apply_adjustments`
   - TLC or no-TLC branch -> `hr1`, `lr1`, `fr1`
   - `create_hr2`

5. **Optional diagnostics and batch bookkeeping**
   - generate overview and section plots
   - refresh `list_tags.csv` and `list_deployments.csv`

## What this means operationally

You can now use the package to:

- check whether a deployment is known
- inspect where raw data and outputs are expected
- stage JSON metadata under `data/data_raw/config_files/`
- place raw ODV and HR files directly under the canonical `data/data_raw/` tree
- process deployments without any external engine selection
- resolve HR raw filenames strictly from `data/catalog/list_deployment_hr.csv`
- rely on batch reruns to prune stale successful entries when canonical outputs have been removed

## HR filename resolution

For FR processing, the raw HR text filename is determined only from the catalog row in `list_deployment_hr.csv`:

`data/data_raw/raw_smru_hr_data/<year>/<prefix><instr_id>_ctd.txt`

The code does not fall back to platform JSON metadata or filename guessing. Packaged sample copies of `list_deployment.csv` and `list_deployment_hr.csv` live under `src/meop_process/resources/catalog/` for tests and examples, but operational runs still read and write the authoritative catalog under `data/catalog/`.

The code also does not fall back to removed root-level runtime locations. Operational runs now use only the canonical `data/` layout.

The remaining deferred work is scientific and analyst-facing, not orchestration-related.
