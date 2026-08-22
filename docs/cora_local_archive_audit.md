# Local CORA tile archive audit

Audit date: 2026-08-20

Configured archive:
`/media/disk2/roquet/REF_DATASETS/CORA_data/CORA_ncfiles`

## Conclusion

The configured directory is a derived, tiled mixture of CORA profiler (`PF`) and CTD
candidate (`CT`) files, but its exact upstream CORA release cannot be established from
the files or this repository. It should not be described as an Argo-only archive or as
a pinned version of the official CORA product.

The files have no global attributes identifying a CORA product code, release, DOI,
creation procedure or source download. No script or manifest that created the tiles is
present in `MEOP_process`. A historical splitter and notebook do exist beside the
archive in `../REF_DATASETS/CORA_data/`; they show an intent to extract both `PF.nc`
and `CT.nc` source files, accept level QC 1 and interpolate to 10--1000 dbar. They do
not identify the download release or prove that the installed files were made by that
exact code. All files share a modification date of 2019-08-12 and the decoded
observations span 2004-01-01 through 2017-12-31. The safest description is:

> Local CORA-derived PF/CT tile archive, generated no later than August 2019 from an
> unidentified upstream release, containing observations through 2017 and lacking a
> reliable retained per-profile source class.

Before this archive is used to justify operational calibration coefficients, it should
be replaced by a build from a named, dated upstream product with source and preparation
manifests. The repository now provides that workflow; see
`docs/cora_download_and_preparation.md`.

## Inventory

| Item | Audit result |
|---|---:|
| NetCDF files | 576 |
| Parsed geographic cells | 576 |
| Duplicate cells | 0 |
| Longitude cells | 36, from 180 W to 170 E |
| Latitude cells | 16, from 80 S to 70 N |
| Populated tiles | 463 |
| Empty header-only placeholders | 113 |
| Profiles in populated tiles | 1,793,760 |
| Archive size | 9,741,037,192 bytes (9.1 GiB on disk) |
| File modification date | 2019-08-12 for all 576 files |
| Decoded observation interval | 2004-01-01 to 2017-12-31 |
| Vertical grid | 100 levels, 10--1000 dbar, identical in populated tiles |
| Profile parameter layouts | 433 tiles with `N_PARAM=3`; 30 with `N_PARAM=4` |
| File-level global attributes | None |
| `DATA_TYPE` | `Argo profile` in all populated tiles, but not reliable after the merge |
| `PROBE_TYPE` | Missing from all 463 populated tiles |
| `WMO_INST_TYPE=830` | 32,228 profiles (CTD instrument code, not a unique platform class) |
| `REFERENCE_DATE_TIME` | `19500101000000` in all populated tiles |

Filename widths confirm the convention that caused the former loader defect: 304 files
use two longitude digits below 100 degrees (for example `lon40W` or `lon00E`), while
272 use three digits because the longitude is at least 100 degrees. The previous loader
always requested three digits.

The scalar `DATA_TYPE` value is contradicted by the recovered splitter and the
per-profile instrument codes. In particular, 32,228 profiles have `WMO_INST_TYPE=830`.
The current CORA manual warns that code 830 can also occur on animal-mounted or moored
CTDs, so it cannot reconstruct the lost platform/source class exactly. Likewise, it is
not safe to label every non-830 profile as Argo. The new preparation workflow uses
`PROBE_TYPE` before merging and retains it in every generated tile.

The filename/size/modification-time manifest fingerprint is:

```text
SHA-256 4f55dcc8ad434ded8c2dd809c157433c43bb8ba9c35cdf229b1b9e7fe6f73e0c
```

This is a lightweight inventory fingerprint, not a content hash of the 9.1 GiB archive.

## Data-mode composition

The populated tiles contain the following `DATA_MODE` values:

| Mode | Profiles | Fraction |
|---|---:|---:|
| Delayed mode (`D`) | 1,061,237 | 59.16% |
| Adjusted real time (`A`) | 496,727 | 27.69% |
| Real time (`R`) | 235,796 | 13.15% |

Consequently, an estimator must define which modes and which raw/adjusted reference
variables are admissible. The archive is not a delayed-mode-only subset.

## Time encoding finding

The 463 populated files contain 453 distinct `JULD` unit origins, generally expressed
as `seconds since <tile-specific date>`. The `REFERENCE_DATE_TIME` variable still says
1950-01-01. This is consistent with a derived xarray-style rewrite rather than untouched
source CORA files.

The current plotting code does not use reference time, so this does not alter its present
figures. It is critical for the proposed space--time matchup work: raw numerical `JULD`
values from different tiles must not be concatenated and compared. Each tile must first
be decoded using its own units/calendar and normalised to UTC.

## Software changes made from this audit

- `reference/cora.py` now inventories and parses the filenames actually present.
- Both `lon40W` and `lon040W` map to the same `(-40, latitude)` cell.
- Duplicate spellings are selected deterministically and logged, never double counted.
- Tests now reproduce the installed archive's unpadded convention.
- A real-archive check for 40--30 W, 70--60 S now loads 1,606 profiles; the old loader
  found no files for that box.
- Calibration batch state now uses method version 4, so incompatible success markers
  are not reused after the target-selection and indexed-loader fixes.
- Invalid/null-island and minor disconnected target components are removed before
  reference selection; impossible targets stop before any CORA tile is opened.
- Reference selection uses the union of local buffered track cells with circular
  anti-meridian handling, rather than a rectangular track envelope.
- Tile coordinates are read first and TEMP/PSAL are materialised only for selected
  profile indices, using the existing 128-profile NetCDF chunks.
- Success now requires at least 10 target profiles and actual T/S support spanning
  200--600 dbar, including at least three levels in the 400--600 dbar band.
- Expected data outcomes are reported separately as `no_reference`,
  `insufficient_target`, `insufficient_reference` and `invalid_target`; only unexpected exceptions are
  `failed`.

## Provenance information still missing from the old archive

1. The official CORA product identifier and release used as input.
2. The acquisition/download date and original filenames.
3. Confirmation that the adjacent historical program, and which exact parameters, made
   these installed files.
4. Whether duplicate removal or additional QC was applied before tiling.
5. Whether any MEOP profiles can occur in the source archive and how to identify them.
6. A content manifest or checksums tied to the source and generated products.

For comparison, the current official Copernicus multi-year discrete CORA product is
[`INSITU_GLO_PHY_TS_DISCRETE_MY_013_001`](https://data.marine.copernicus.eu/product/INSITU_GLO_PHY_TS_DISCRETE_MY_013_001/services);
the Marine Data Store currently exposes dataset version `202511`. That fact does not
identify this 2019 local archive as any particular release of that product. The
[official manual](https://documentation.marine.copernicus.eu/PUM/CMEMS-INS-PUM-013-001.pdf)
identifies CTD and profiler classes with `PROBE_TYPE=2` and `PROBE_TYPE=5` respectively,
and explicitly warns that file class or WMO instrument code alone is not sufficient for
every CTD.

## Audit method

The audit parsed every filename with the CORA longitude/latitude pattern and opened all
576 NetCDF headers. For populated files it counted `N_PROF`, compared variable schemas
and `PRES_GRID`, decoded `JULD` with each file's own units/calendar, and counted
`DATA_MODE`. A follow-up pass counted `PROBE_TYPE` and `WMO_INST_TYPE` and reviewed the
historical code adjacent to the archive. `ncdump -h` confirmed the absence of global
provenance attributes. No source download manifest was found.
