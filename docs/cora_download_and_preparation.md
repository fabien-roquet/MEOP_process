# Downloading and preparing a current CTD/Argo CORA reference archive

## Short answer

Yes. The installed local tiles stop on 2017-12-31, so they omit every profile from
2018 onward and do not include later reprocessing of older observations. The current
Copernicus Marine file browser exposes CORA dataset version `202511`; CORA is updated
more frequently for Argo than for other platforms. Do not infer the exact last profile
date from the release label: the preparation manifest and a data-level audit should
record the actual minimum and maximum times after downloading.

The old local tiles cannot be separated reliably into CTD and Argo profiles. The
historical script extracted both `PF` and `CT` files, but discarded `PROBE_TYPE` while
merging. All populated tiles misleadingly retain `DATA_TYPE="Argo profile"`, while
32,228 profiles have WMO instrument code 830 (CTD). The official manual also warns that
instrument code 830 is not unique to ship CTDs. It is therefore unsafe either to treat
the archive as Argo-only or to classify every profile from `WMO_INST_TYPE` alone.

The replacement workflow makes the distinction before merging:

- CTD: raw CORA `CT`, legacy `OC`, and legacy `TE` candidate files, followed by the
  per-profile condition `PROBE_TYPE=2`;
- Argo by default: EasyCORA `PF` files, which the official manual specifically describes
  as profiles received from Argo DACs;
- optional full-level alternative: raw CORA `PF` plus `PROBE_TYPE=5`. CORA calls code 5
  “profilers,” so this is not as strict a statement of Argo-program membership.

The default EasyCORA Argo branch contains CORA's best parameter estimate and has already
been vertically subsampled. This is appropriate for the 10 dbar MEOP comparison tiles,
but the source choice is written into every profile and manifest. Use the raw-PF option
if preserving original levels is more important than strict Argo provenance.

Official references:

- [CORA product and current data access](https://data.marine.copernicus.eu/product/INSITU_GLO_PHY_TS_DISCRETE_MY_013_001/services)
- [CORA Product User Manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-INS-PUM-013-001.pdf)
- [Copernicus Marine `get` documentation](https://help.marine.copernicus.eu/en/articles/8286883-copernicus-marine-toolbox-api-get-original-files)

## 1. Install and authenticate

From the repository root, use a Python environment with enough disk space:

```bash
python -m pip install -e '.[dev,cora]'
copernicusmarine login
```

Credentials are stored by the Copernicus Marine Toolbox. The repository scripts do not
accept or write passwords.

## 2. Plan before downloading

Use a new source directory for each upstream release:

```bash
meop-cora-download plan \
  --output-dir /path/to/references/cora_ctd_argo_source
```

This does not download profile NetCDF files. It creates:

```text
cora_ctd_argo_source/
  download_manifest.json
  manifests/
    <dataset>_latest_remote_inventory.csv
    selected_ctd.csv
    selected_ctd.txt
    selected_argo.csv
    selected_argo.txt
```

The command prints the selected file counts, reported byte totals and discovered
dataset version. Inspect `download_manifest.json` and confirm the storage requirement.
The exact lists and their SHA-256 hashes make the acquisition reproducible.

The default selects the complete time series. This is recommended because newer CORA
releases can reprocess old profiles. To make a deliberately incomplete test download:

```bash
meop-cora-download plan \
  --output-dir /path/to/references/cora_test_source \
  --start-year 2024 --end-year 2024
```

Do not append such a recent-only selection to the unidentified old tiles and call the
result one archive: duplicate/QC/value policies would differ.

To use original-level raw profiler files instead of the strict EasyCORA Argo stream,
add `--argo-source raw-cora`. If a release must be pinned in advance, add for example
`--dataset-version 202511`; otherwise the latest version discovered from the remote
paths is pinned automatically in the saved manifest.

## 3. Download the exact planned files

Run the same selection with `download`:

```bash
meop-cora-download download \
  --output-dir /path/to/references/cora_ctd_argo_source
```

Existing files are skipped. Sources are kept separate under:

```text
source/cora_ctd_candidates/
source/easycora_argo/
```

For a later CORA release, prefer a new source root. If intentionally refreshing a saved
`latest` inventory in place, add `--refresh-manifest`, inspect the changed plan, and only
then download.

## 4. Prepare MEOP-compatible tiles

Build into a new output directory rather than replacing the installed 2019 archive:

```bash
meop-cora-prepare \
  --source-root /path/to/references/cora_ctd_argo_source \
  --output-dir /path/to/references/CORA_ctd_argo_202511_tiles
```

The preparation command:

1. enforces source class and `PROBE_TYPE` before merging;
2. rejects bad position/time QC when those flags exist;
3. uses adjusted values for a whole parameter profile when available, otherwise raw;
4. retains only level QC 1 (good) and 2 (probably good);
5. requires both usable temperature and salinity by default;
6. interpolates without extrapolation to 10--1000 dbar every 10 dbar;
7. normalises longitude and time, writing one common `JULD` origin;
8. writes 10-degree tiles accepted by `reference/cora.py`;
9. preserves profile kind, probe type, identifiers, source file/class/dataset, and the
   raw/adjusted branch used for temperature, salinity, and pressure.

Use `--allow-single-variable` only when temperature-only or salinity-only profiles are
needed. The command stops on malformed files by default; `--skip-bad-files` is intended
for diagnosis and records every skip in the manifest. Existing tiles are never replaced
unless the exact output directory is supplied with `--overwrite`.

The result contains `PREPARATION_MANIFEST.json`. Review at least:

- source inventory fingerprint and source file count;
- `profiles_written_argo` and `profiles_written_ctd`;
- type-matching profiles rejected for time/position or insufficient T/S;
- missing QC variables and skipped-file errors;
- pressure grid and selection/value/QC policies.

## 5. Configure and make a small validation run

Point `references.cora_dir` in `configs.json` to the new tile directory. Keep the old
path available until the new archive has passed review.

```json
{
  "defaults": {
    "references": {
      "cora_dir": "/path/to/references/CORA_ctd_argo_202511_tiles"
    }
  }
}
```

Then run a small test:

```bash
meop-compare --plot1 ct88-225-12
meop-compare-batch --dry-run --deployment ct88
pytest -q tests/test_cora_download_prepare.py tests/test_cora_calibration.py
```

Before a scientific batch, audit actual coverage dates and counts by reference kind,
check several raw source profiles against generated tiles, test a dateline region, and
confirm that MEOP/animal-mounted profiles can be excluded using retained identifiers.
The new archive is suitable infrastructure for matchups; it does not by itself validate
an adjustment estimator.
