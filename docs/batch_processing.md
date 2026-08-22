# Batch processing and resumable reruns

The toolbox includes a Python-native batch runner for processing many deployments in sequence.

## Goals

The runner is designed to:

- continue past deployment-level failures;
- avoid redoing deployments that already completed successfully when canonical outputs still exist;
- keep one readable log per deployment;
- produce a compact Markdown and CSV summary for each batch run;
- refresh `list_tags.csv` and `list_deployments.csv` at the end of the run.

## Main commands

Editable install:

```bash
meop-process-batch
```

Main CLI entrypoint:

```bash
meop-process --run-all-deployments
```

Repository checkout:

```bash
python scripts/run_all_deployments.py
```

## Common options

Diagnostics are enabled by default for successful deployments (`defaults.batch.diagnostics: true`):

This writes the standard per-tag figures under `data/plots_by_tags/`, a deployment recap figure under `data/plots_by_deployments/` for each processed deployment, and cross-deployment overview summaries under `data/plots_overview/` when diagnostics are run across multiple deployments.

```bash
meop-process --run-all-deployments
meop-process-batch
```

Disable diagnostics explicitly for one run:

```bash
meop-process --run-all-deployments --no-diagnostics
meop-process-batch --no-diagnostics
```

Successful production runs keep only the configured final products (`hr1`, `lr1`, and `hr2` by default) and prune rebuildable intermediates such as `lr0`, `hr0`, `fr0`, and `fr1`.
Keep those stage products for calibration/debug runs with:

```bash
meop-process --run-all-deployments --keep-intermediate-products
meop-process-batch --keep-intermediate-products
```

Restrict diagnostics to one layer:

```bash
meop-process --run-all-deployments --diagnostics --diagnostics-part overview
meop-process --run-all-deployments --diagnostics --diagnostics-part deployment
meop-process --run-all-deployments --diagnostics --diagnostics-part tag
```

Send the final batch summary by email:

```bash
meop-process --run-all-deployments --diagnostics --notify-email ops@example.org
```

The batch summary email can also be enabled in root `configs.json` under `defaults.notifications.email` or `configs.<machine>.notifications.email`.

```bash
python scripts/run_all_deployments.py --diagnostics
```

Re-run only deployments whose latest status is failed:

```bash
meop-process --run-all-deployments --force-failed
```

```bash
python scripts/run_all_deployments.py --force-failed
```

Force a full rerun even for deployments that already completed successfully:

```bash
meop-process --run-all-deployments --force
```

## CORA-based T/S calibration plots

After processing, CORA-based T/S calibration plots can be generated for a single tag using `meop-compare`:

```bash
meop-compare --plot1 ct88-225-12
```

This command does not compare two output trees. Instead it:

1. Reads the processed profile file selected for the tag.
2. Removes invalid/null-island positions and minor spatial components disconnected
   from a supported deployment component.
3. Rejects targets that cannot meet the profile/depth gate before opening CORA.
4. Loads only the union of CORA cells intersecting local 5° windows along the accepted
   track, with circular longitude handling at the anti-meridian.
5. Reads tile coordinates first and materialises TEMP/PSAL only for selected profile rows.
6. Writes one PNG per 200-profile chunk to `data/plots_by_tags/<deployment>/`.

The tile loader inventories actual filenames and supports both `lon40W`/`lon00E` and
`lon040W`/`lon000E` conventions. Batch summaries classify `no_reference`,
`insufficient_target`, `insufficient_reference`, and `invalid_target` separately from
unexpected `failed` outcomes. Calibration method version 4 invalidates incompatible
success state. A tag is
successful only when both temperature and salinity diagnostics contain at least 10
target profiles, at least three valid levels in the 400--600 dbar band, finite band and
linear coefficient estimates, and observed comparison coverage from 200 to 600 dbar.

Each figure has two panels:
- **Left** — T/S diagram with CORA background profiles (grey), other tags in the same deployment (blue), and the target tag coloured by time (viridis).
- **Right** — PSAL anomaly versus pressure relative to the per-level CORA median.

### Required configuration

`cora_dir` must be set in root `configs.json`:

```json
{
  "defaults": {
    "references": {
      "cora_dir": "/path/to/CORA_ncfiles"
    }
  },
  "configs": {
    "my_machine": {
      "processdir": "/path/to/MEOP_process",
      "datadir": "data",
      "public": "public"
    }
  }
}
```

If `cora_dir` is not configured, the command exits with an informative error message.

A custom config file can also be supplied directly:

```bash
meop-compare --plot1 ct88-225-12 --config /path/to/configs.json
```


```bash
python scripts/run_all_deployments.py --force
```

Restrict a batch to one deployment:

```bash
meop-process --run-all-deployments --deployment ct96
```

```bash
python scripts/run_all_deployments.py --deployment ct96
```

Use the no-TLC branch:

```bash
meop-process --run-all-deployments --notlc
```

```bash
python scripts/run_all_deployments.py --notlc
```

Run several deployments in parallel and mirror deployment logs to the terminal:

```bash
meop-process --run-all-deployments --jobs 8 --verbose
```

## State, logs, and reports

By default, batch artifacts are stored under `data/batch/`:

- `latest/deployment_status.json`: latest persistent state keyed by deployment code
- `runs/<timestamp>/logs/<deployment>.log`: per-deployment logs
- `runs/<timestamp>/summary.md`: human-readable summary
- `runs/<timestamp>/summary.csv`: machine-readable summary

The state file is what makes the runner resumable.
A deployment marked `success` is skipped on the next run only if its canonical outputs still exist under `data/data_prof/`, unless `--force` is used.
A batch run reconciles `latest/deployment_status.json` against the filesystem at startup and prunes stale successful entries whose outputs have been deleted.
Deployments marked `failed` are attempted again on the next run.

## Metadata summary refresh

At the end of the batch run, `list_tags.csv` and `list_deployments.csv` are refreshed.
The update is incremental:

- processed deployments are refreshed;
- deployments whose output tag inventory changed are refreshed;
- unchanged deployments are preserved as-is.

This keeps summary refresh fast on reruns.
