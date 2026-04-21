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

Run diagnostics for successful deployments:

This writes the standard per-tag figures under `data/plots_by_tags/`, a deployment recap figure under `data/plots_by_deployments/` for each processed deployment, and cross-deployment overview summaries under `data/plots_overview/` when diagnostics are run across multiple deployments.

```bash
meop-process --run-all-deployments --diagnostics
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

The batch summary email can also be enabled in `data/configs.json` under `defaults.notifications.email` or `configs.<machine>.notifications.email`.

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
