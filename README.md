# MEOP_process
Scripts used to process MEOP data (meop.net)

Since 2004, several hundred thousands profiles of temperature and salinity have been 
collected by instrumented animals. The use of elephant seals has been particularly 
effective to sample the Southern Ocean and the North Pacific. These hydrographic data 
have been assembled in a quality-controlled database, the MEOP-CTD database, that can 
be accessed through this website.
For more information, visit the website meop.net
For any questions, contact info@meop.net or fabien.roquet@gmail.com


## THE MEOP-CTD DATABASE
README: the MEOP-CTD database (owner: Fabien Roquet and the MEOP consortium)
Release Date: 11/11/2017
Version name: MEOP-CTD_2017-11-11
https://opendatacommons.org/licenses/odbl/


## DATA FORMATS

Data are provided in three different formats. 
* _DATA_ncARGO:_ For a thourough scientific use of the data, or for oceanographic data centers, it is advised to use the 
marine mammal netCDF format (files in DATA_ncARGO) as it serves as the reference. This format can be 
easily read in Ocean Data View, using the Import/ARGO profiles/Float profiles menu, or using your 
favorite data processing software (e.g. Python, R, IDL). 
* _DATA_ncARGO_interp:_ For ease of use, the DATA_ncARGO_interp provides the same data as in DATA_ncARGO, except it has
been interpolated on a regular vertical grid (1dbar spacing).
* _DATA_csv_interp:_ A csv format (ASCII) is also provided (files in DATA_csv_interp) which can be opened with Excel
or any text editor. Here, only data flagged as good are included, and are given on a regular 
vertical grid (1dbar spacing).



## DATA INFORMATION

The data that is publicly available is shown in the figure map_global_public.png. More data 
is available upon request, as part of the MEOP-CTD database. See map_global_private.png for
the distribution of private data.
Important metadata and statistics are listed in the info_*.csv files:
* _info_total.csv_ gives global statistics about the MEOP-CTD database
* _info_groups.csv_ gives statistics by national groups (see [MEOP groups](meop.net/groups/) for information)
* _info_deployments.csv_ gives statistics by deployment
* _info_tags.csv_ gives information and statistics by individual tag.
For each deployments, distribution maps are available in the MAPS directory, and a pdf
document with basic plots of CTD data (TS plots, time sections) is provided in the
directory PDF.



## HOW TO CITE

If you use this dataset for a publication, please add the following sentence 
in the Acknowledgement part:
"The marine mammal data were collected and made freely available by the International MEOP 
Consortium and the national programs that contribute to it (http://www.meop.net)."

Also consider citing the following papers when you use the MEOP-CTD dataset
for oceanographic applications:
- Treasure, A. M., Roquet, F., Ansorge, I. J., Bester, M. N., Horst Bornemann, L. B., Charrassin, J.-B., Chevallier, D., Costa, D. P., Fedak, M. A., Guinet, C., Hammill, M. O., Harcourt, R. G., Hindell, M. A., Kovacs, K. M., Lea, M.-A., Lovell, P., Lowther, A. D., Lydersen, C., McIntyre, T., McMahon, C. R., Muelbert, M. M. C., Nicholls, K., Picard, B., Reverdin, G., Trites, A. W., Williams, G. D., and de Bruyn, P. J. N., 2017. Marine Mammals Exploring the Oceans Pole to Pole: A Review of the MEOP Consortium. Oceanography, 30(2):132–138, doi: 10.5670/oceanog.2017.234
- Roquet F., Williams G., Hindell M. A., Harcourt R., McMahon C. R., Guinet C., Charrassin 
J.-B., Reverdin G., Boehme L., Lovell P. and Fedak M. A., 2014. A Southern Indian Ocean 
database of hydrographic profiles obtained with instrumented elephant seals. Nature 
Scientific Data, 1:140028, doi: 10.1038/sdata.2014.28

### Important technical papers : 
* A thorough description of the CTD-SRDL technology can be found in : 
  * Boehme L., Lovell P., Biuw M., Roquet F., Nicholson J., Thorpe S. E., Meredith M. P., and Fedak M., 2009. Technical Note: Animal-borne CTD-Satellite Relay Data Loggers for real-time oceanographic data collection. Ocean Sci., 5:685-695. doi: 10.5194/os-5-685-2009
* The delayed-mode general methodology and estimated accuracy of CTD-SRDL hydrographic data 
are presented in :
  * Roquet F., Charrassin J.-B., Marchand S., Boehme L., Fedak M., Reverdin G., and Guinet C., 2011. Validating hydrographic data obtained from seal-borne satellite-relayed data loggers. J. Atmos. Oceanic Technol., 28:787-801. doi: 10.1175/2010JTECHO801.1
* The density inversion removal algorithm is described in :
  * Barker, P. M. and McDougall, T. J., 2017. Stabilizing Hydrographic Profiles with Minimal Change to the Water Masses. J. Atmos. Oceanic Technol., 34:1935-1945. doi: 10.1175/JTECH-D-16-0111.1
* The thermal cell effect correction is described in :
  * Mensah, V., Roquet, F., Siegelman-Charbit, L., Picard, B., Pauthenet, E., Guinet, C., 2018.  A correction methodology for the thermal mass induced-errors of CTD tags mounted on marine mammals. J. Atmos. Oceanic Technol., 35:1237–1252. doi: 10.1175/JTECH-D-17-0141.1


### National specificities :
- For Australian data: 
Any users of IMOS data are required to clearly acknowledge the source of the material 
derived from IMOS in the format: 
"Data was sourced from the Integrated Marine Observing System (IMOS) - IMOS is a national 
collaborative research infrastructure, supported by Australian Government.” IMOS data is 
licensed under a Creative Commons Attribution (CCBY) License, 
(http://creativecommons.org.au/)."
- For German and South African data: 
Primary data are also made available through PANGAEA. Please cite :
doi10.1594/PANGAEA.150008 for data related to Marion Island (Southern Ocean Indian Sector)
doi10.1594/PANGAEA.150009 for data related to King George Island (Southern Ocean Atlantic Sector)
doi10.1594/PANGAEA.150010 for data related to Atka Bay, Drescher Inlet, Filchner Trough (Southern Ocean Atlantic Sector)

## CURRENT SOFTWARE STATUS

The repository now contains a pure-Python package under `src/meop_process/`.
The current minimal functional pipeline covers:

- deployment and tag discovery from catalog CSV files and JSON metadata;
- raw ODV import and profile indexing;
- `lr0`, QC/filtering, `hr0`, `hr1`, `lr1`, `fr0`, `fr1`, and `hr2`;
- delayed-mode `apply_adjustments`;
- standard diagnostics figures;
- batch processing over multiple deployments with resumable state, readable reports, and per-deployment logs.

The package now has a single Python execution path and no longer requires root-level table or catalog mirrors.

## INSTALLATION

A typical editable install is:

```bash
python -m pip install -e .
```

The current Python runtime expects at least:

- `numpy`
- `pandas`
- `xarray`
- `scipy`
- `h5netcdf`
- `h5py`
- `matplotlib`
- `gsw`

If cartographic map backgrounds are desired in diagnostics, install `cartopy` as an additional optional dependency.

## RUNTIME DATA LAYOUT

The cleaned package expects data in explicit runtime locations:

- packaged/default tables: `src/meop_process/resources/tables/`, synchronized into `data/tables/`
- operator-managed catalog tables: `data/catalog/`
- deployment/platform JSON files: `data/config_files/`
- raw low-resolution ODV files: `data/raw_smru_data_odv/`
- raw high-resolution text files: `data/raw_smru_hr_data/<year>/<instr_id>_ctd.txt`
- references and comparison datasets: `references/`
- processed outputs: `final_dataset_prof/`
- diagnostics: `plots/`
- batch logs and resumable state: `data/batch/`

## RUNNING ONE DEPLOYMENT

From an editable install:

```bash
meop-process --deployment ct88 --process_data --diagnostics
```

From the repository checkout:

```bash
python python/meop_process.py --deployment ct88 --process_data --diagnostics
```

## BATCH RERUN OVER ALL DEPLOYMENTS

A resumable batch runner is available.
It continues past errors, writes one log per deployment, generates a readable Markdown summary plus a CSV report, and does not redo successful deployments unless forced.

Installed entry point:

```bash
meop-process-batch
```

Repository wrapper script:

```bash
python scripts/run_all_deployments.py
```

Useful options:

```bash
python scripts/run_all_deployments.py --diagnostics
python scripts/run_all_deployments.py --force-failed
python scripts/run_all_deployments.py --force
python scripts/run_all_deployments.py --deployment ct96
python scripts/run_all_deployments.py --notlc
```

Batch state and reports are stored under `data/batch/` by default:

- `data/batch/latest/deployment_status.json`: latest persistent per-deployment state
- `data/batch/runs/<timestamp>/logs/`: one log file per deployment
- `data/batch/runs/<timestamp>/summary.md`: human-readable run report
- `data/batch/runs/<timestamp>/summary.csv`: machine-readable run table

## METADATA SUMMARY TABLES

At the end of a batch run, the package refreshes `list_tags.csv` and `list_deployments.csv`.
This update is incremental:

- deployments processed in the current run are refreshed;
- deployments whose tag inventory changed are refreshed;
- unchanged deployments are preserved without reopening their netCDF files.

The output directory is resolved automatically:

- if an existing `list_tags.csv` / `list_deployments.csv` already exists under the configured public root, it is updated in place;
- otherwise the files are written under `public/<version>/`.

You can also refresh those summary CSVs without reprocessing deployments:

```bash
meop-process --refresh-metadata-summaries
```
