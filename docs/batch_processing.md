# Batch processing and resumable reruns

The toolbox includes a Python-native batch runner for processing many deployments in sequence.

## Goals

The runner is designed to:

- continue past deployment-level failures;
- avoid redoing deployments that already completed successfully;
- keep one readable log per deployment;
- produce a compact Markdown and CSV summary for each batch run;
- refresh `list_tags.csv` and `list_deployments.csv` at the end of the run.

## Main commands

Editable install:

```bash
meop-process-batch
```

Repository checkout:

```bash
python scripts/run_all_deployments.py
```

## Common options

Run diagnostics for successful deployments:

```bash
python scripts/run_all_deployments.py --diagnostics
```

Re-run only deployments whose latest status is failed:

```bash
python scripts/run_all_deployments.py --force-failed
```

Force a full rerun even for deployments that already completed successfully:

```bash
python scripts/run_all_deployments.py --force
```

Restrict a batch to one deployment:

```bash
python scripts/run_all_deployments.py --deployment ct96
```

Use the no-TLC branch:

```bash
python scripts/run_all_deployments.py --notlc
```

## State, logs, and reports

By default, batch artifacts are stored under `data/batch/`:

- `latest/deployment_status.json`: latest persistent state keyed by deployment code
- `runs/<timestamp>/logs/<deployment>.log`: per-deployment logs
- `runs/<timestamp>/summary.md`: human-readable summary
- `runs/<timestamp>/summary.csv`: machine-readable summary

The state file is what makes the runner resumable.
A deployment marked `success` is skipped on the next run unless `--force` is used.
Deployments marked `failed` are attempted again on the next run.

## Metadata summary refresh

At the end of the batch run, `list_tags.csv` and `list_deployments.csv` are refreshed.
The update is incremental:

- processed deployments are refreshed;
- deployments whose output tag inventory changed are refreshed;
- unchanged deployments are preserved as-is.

This keeps summary refresh fast on reruns.
