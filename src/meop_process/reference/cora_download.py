"""Plan and download a traceable CTD/Argo subset of CORA.

The downloader deliberately creates a complete remote inventory first.  It then
selects only the global files needed for ship CTDs and Argo profiles and saves the
exact selection before asking the Copernicus Marine Toolbox to download anything.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Sequence


CORA_DATASET_ID = "cmems_obs-ins_glo_phy-temp-sal_my_cora_irr"
EASYCORA_DATASET_ID = "cmems_obs-ins_glo_phy-temp-sal_my_easycora_irr"
PRODUCT_ID = "INSITU_GLO_PHY_TS_DISCRETE_MY_013_001"

_PROFILE_NAME = re.compile(
    r"_(?P<date>\d{8})_PR_(?P<file_class>[A-Z0-9]{2})\.nc$",
    re.IGNORECASE,
)
_DATASET_VERSION = re.compile(r"_(?P<version>20\d{4})(?:/|$)")


@dataclass(frozen=True)
class RemoteFile:
    """One row from a Copernicus Marine ``get`` CSV inventory."""

    filename: str
    size: int | None = None
    last_modified_datetime: str = ""
    etag: str = ""


@dataclass(frozen=True)
class SelectionSpec:
    """Rules for one independently downloadable source selection."""

    name: str
    dataset_id: str
    file_classes: tuple[str, ...]
    destination: str
    description: str


def selection_specs(argo_source: str) -> tuple[SelectionSpec, ...]:
    """Return the source selections needed for the requested Argo policy."""
    ctd = SelectionSpec(
        name="ctd",
        dataset_id=CORA_DATASET_ID,
        file_classes=("CT", "OC", "TE"),
        destination="cora_ctd_candidates",
        description=(
            "Raw CORA profile files that can contain CTDs; preparation keeps only "
            "profiles with PROBE_TYPE=2."
        ),
    )
    if argo_source == "easycora":
        argo = SelectionSpec(
            name="argo",
            dataset_id=EASYCORA_DATASET_ID,
            file_classes=("PF",),
            destination="easycora_argo",
            description=(
                "EasyCORA PF files, documented as Argo profiles received from the "
                "Argo DACs; values are best-estimate and vertically subsampled."
            ),
        )
    elif argo_source == "raw-cora":
        argo = SelectionSpec(
            name="argo",
            dataset_id=CORA_DATASET_ID,
            file_classes=("PF",),
            destination="cora_argo",
            description=(
                "Raw CORA PF files; preparation keeps PROBE_TYPE=5. CORA calls this "
                "the profiler class, so it is less strict than EasyCORA PF provenance."
            ),
        )
    else:  # pragma: no cover - argparse prevents this for CLI callers
        raise ValueError(f"Unsupported Argo source: {argo_source}")
    return ctd, argo


def read_remote_inventory(path: Path | str) -> list[RemoteFile]:
    """Read the CSV produced by ``copernicusmarine.get(create_file_list=...)``."""
    inventory_path = Path(path)
    with inventory_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Remote inventory has no header: {inventory_path}")
        path_column = next(
            (name for name in ("filename", "file_path", "path") if name in reader.fieldnames),
            None,
        )
        if path_column is None:
            raise ValueError(
                f"Remote inventory {inventory_path} has no filename column; "
                f"found {reader.fieldnames}"
            )
        rows: list[RemoteFile] = []
        for raw in reader:
            filename = str(raw.get(path_column, "")).strip()
            if not filename:
                continue
            size_raw = str(raw.get("size", "")).strip()
            try:
                size = int(size_raw) if size_raw else None
            except ValueError:
                size = None
            rows.append(
                RemoteFile(
                    filename=filename,
                    size=size,
                    last_modified_datetime=str(raw.get("last_modified_datetime", "")).strip(),
                    etag=str(raw.get("etag", "")).strip(),
                )
            )
    return rows


def _profile_details(filename: str) -> tuple[int, str] | None:
    match = _PROFILE_NAME.search(Path(filename).name)
    if match is None:
        return None
    return int(match.group("date")[:4]), match.group("file_class").upper()


def _has_global_component(filename: str) -> bool:
    return any(component.casefold() == "global" for component in filename.split("/"))


def select_remote_files(
    rows: Iterable[RemoteFile],
    *,
    file_classes: Iterable[str],
    start_year: int | None = None,
    end_year: int | None = None,
) -> list[RemoteFile]:
    """Select global daily profile files by class and inclusive year range.

    Requiring the explicit ``Global`` path component prevents downloading duplicate
    regional extractions of the same product.
    """
    allowed = {item.upper() for item in file_classes}
    candidates: list[RemoteFile] = []
    for row in rows:
        details = _profile_details(row.filename)
        if details is None:
            continue
        year, file_class = details
        if file_class not in allowed:
            continue
        if start_year is not None and year < start_year:
            continue
        if end_year is not None and year > end_year:
            continue
        candidates.append(row)

    selected = [row for row in candidates if _has_global_component(row.filename)]
    if candidates and not selected:
        raise ValueError(
            "Matching profile files were found, but none had an explicit Global path "
            "component. Refusing to mix regional extractions; inspect the saved remote "
            "inventory and update the selector for the new upstream layout."
        )
    return sorted(selected, key=lambda item: item.filename)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_selection(rows: Sequence[RemoteFile], txt_path: Path, csv_path: Path) -> None:
    txt_path.write_text("".join(f"{row.filename}\n" for row in rows), encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("filename", "size", "last_modified_datetime", "etag"),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def _versions(rows: Iterable[RemoteFile]) -> list[str]:
    versions: set[str] = set()
    for row in rows:
        for match in _DATASET_VERSION.finditer(row.filename):
            versions.add(match.group("version"))
    return sorted(versions)


def _load_toolbox() -> Any:
    try:
        import copernicusmarine  # type: ignore
    except ImportError as exc:  # pragma: no cover - depends on optional install
        raise RuntimeError(
            "The optional Copernicus Marine Toolbox is not installed. Run "
            "`python -m pip install -e '.[cora]'` from this repository, then "
            "`copernicusmarine login`."
        ) from exc
    return copernicusmarine


def _inventory_dataset(
    toolbox: Any,
    *,
    dataset_id: str,
    dataset_version: str | None,
    manifest_dir: Path,
    refresh: bool,
) -> Path:
    version_tag = dataset_version or "latest"
    csv_path = manifest_dir / f"{dataset_id}_{version_tag}_remote_inventory.csv"
    if csv_path.exists() and not refresh:
        return csv_path
    if csv_path.exists():
        csv_path.unlink()
    kwargs: dict[str, Any] = {
        "dataset_id": dataset_id,
        "output_directory": manifest_dir,
        "create_file_list": csv_path.name,
    }
    if dataset_version:
        kwargs["dataset_version"] = dataset_version
    toolbox.get(**kwargs)
    if not csv_path.exists():
        raise RuntimeError(f"Copernicus Marine did not create {csv_path}")
    return csv_path


def create_download_plan(
    output_dir: Path | str,
    *,
    argo_source: str = "easycora",
    dataset_version: str | None = None,
    start_year: int | None = None,
    end_year: int | None = None,
    refresh: bool = False,
    toolbox: Any | None = None,
) -> dict[str, Any]:
    """Inventory the current upstream release and write exact selected file lists."""
    if start_year is not None and end_year is not None and start_year > end_year:
        raise ValueError("start_year must be less than or equal to end_year")
    root = Path(output_dir).expanduser().resolve()
    manifest_dir = root / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    client = toolbox if toolbox is not None else _load_toolbox()

    specs = selection_specs(argo_source)
    inventory_paths: dict[str, Path] = {}
    for dataset_id in sorted({spec.dataset_id for spec in specs}):
        inventory_paths[dataset_id] = _inventory_dataset(
            client,
            dataset_id=dataset_id,
            dataset_version=dataset_version,
            manifest_dir=manifest_dir,
            refresh=refresh,
        )

    plan_selections: list[dict[str, Any]] = []
    for spec in specs:
        inventory_path = inventory_paths[spec.dataset_id]
        selected = select_remote_files(
            read_remote_inventory(inventory_path),
            file_classes=spec.file_classes,
            start_year=start_year,
            end_year=end_year,
        )
        if not selected:
            raise RuntimeError(
                f"No {spec.name} files matched {spec.file_classes} in {inventory_path}"
        )
        basenames = [Path(row.filename).name for row in selected]
        seen_names: set[str] = set()
        duplicate_name_set: set[str] = set()
        for name in basenames:
            if name in seen_names:
                duplicate_name_set.add(name)
            else:
                seen_names.add(name)
        duplicate_names = sorted(duplicate_name_set)
        if duplicate_names:
            raise RuntimeError(
                f"Selection {spec.name} contains duplicate basenames, which cannot be "
                f"downloaded safely with no_directories=True: {duplicate_names[:5]}"
            )
        txt_path = manifest_dir / f"selected_{spec.name}.txt"
        csv_path = manifest_dir / f"selected_{spec.name}.csv"
        _write_selection(selected, txt_path, csv_path)
        versions = _versions(selected)
        if len(versions) > 1:
            raise RuntimeError(
                f"Selection {spec.name} spans multiple dataset versions: {versions}"
            )
        plan_selections.append(
            {
                "name": spec.name,
                "dataset_id": spec.dataset_id,
                "dataset_version": versions[0] if versions else dataset_version,
                "file_classes": list(spec.file_classes),
                "description": spec.description,
                "destination": str(Path("source") / spec.destination),
                "selected_file_list": str(txt_path.relative_to(root)),
                "selected_metadata": str(csv_path.relative_to(root)),
                "selected_file_list_sha256": _sha256(txt_path),
                "file_count": len(selected),
                "known_bytes": sum(row.size for row in selected if row.size is not None),
                "files_without_reported_size": sum(row.size is None for row in selected),
            }
        )

    plan: dict[str, Any] = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "product_id": PRODUCT_ID,
        "argo_source": argo_source,
        "requested_dataset_version": dataset_version,
        "start_year": start_year,
        "end_year": end_year,
        "remote_inventories": {
            dataset_id: {
                "path": str(path.relative_to(root)),
                "sha256": _sha256(path),
            }
            for dataset_id, path in inventory_paths.items()
        },
        "selections": plan_selections,
    }
    plan_path = root / "download_manifest.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return plan


def download_plan(
    output_dir: Path | str,
    plan: dict[str, Any],
    *,
    toolbox: Any | None = None,
) -> None:
    """Download each exact selection in a previously generated plan."""
    root = Path(output_dir).expanduser().resolve()
    client = toolbox if toolbox is not None else _load_toolbox()
    for selection in plan["selections"]:
        file_list_path = root / selection["selected_file_list"]
        if _sha256(file_list_path) != selection["selected_file_list_sha256"]:
            raise RuntimeError(f"Selected file list changed after planning: {file_list_path}")
        destination = root / selection["destination"]
        destination.mkdir(parents=True, exist_ok=True)
        kwargs: dict[str, Any] = {
            "dataset_id": selection["dataset_id"],
            "file_list": file_list_path,
            "output_directory": destination,
            "no_directories": True,
            "skip_existing": True,
        }
        if selection.get("dataset_version"):
            kwargs["dataset_version"] = selection["dataset_version"]
        client.get(**kwargs)
        selected_rows = read_remote_inventory(root / selection["selected_metadata"])
        missing: list[str] = []
        size_mismatches: list[str] = []
        local_bytes = 0
        for row in selected_rows:
            local_path = destination / Path(row.filename).name
            if not local_path.is_file():
                missing.append(local_path.name)
                continue
            size = local_path.stat().st_size
            local_bytes += size
            if row.size is not None and size != row.size:
                size_mismatches.append(
                    f"{local_path.name}: expected {row.size}, found {size}"
                )
        selection["download"] = {
            "checked_at_utc": datetime.now(timezone.utc).isoformat(),
            "local_file_count": len(selected_rows) - len(missing),
            "local_bytes": local_bytes,
            "missing_files": missing,
            "size_mismatches": size_mismatches,
        }
        if missing or size_mismatches:
            (root / "download_manifest.json").write_text(
                json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            raise RuntimeError(
                f"Download verification failed for {selection['name']}: "
                f"{len(missing)} missing, {len(size_mismatches)} size mismatches"
            )
    plan["download_completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    (root / "download_manifest.json").write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _format_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024.0 or unit == "TiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024.0
    return f"{amount:.1f} TiB"  # pragma: no cover


def _print_plan(root: Path, plan: dict[str, Any]) -> None:
    print(f"Plan written to {root / 'download_manifest.json'}")
    for item in plan["selections"]:
        version = item.get("dataset_version") or "latest (not encoded in paths)"
        print(
            f"{item['name']}: {item['file_count']} files, "
            f"{_format_bytes(item['known_bytes'])}, dataset version {version}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plan or download a traceable CTD/Argo-only CORA source archive."
    )
    parser.add_argument("command", choices=("plan", "download"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("references/cora_ctd_argo_source"),
        help="Source archive root (default: references/cora_ctd_argo_source).",
    )
    parser.add_argument(
        "--argo-source",
        choices=("easycora", "raw-cora"),
        default="easycora",
        help=(
            "Use strict Argo-DAC EasyCORA PF files (default), or full-level raw CORA "
            "PF/profiler files."
        ),
    )
    parser.add_argument(
        "--dataset-version",
        default=None,
        help="Optional YYYYMM release pin; otherwise inventory the latest release and pin it in the manifest.",
    )
    parser.add_argument("--start-year", type=int, default=None, help="First observation year (inclusive).")
    parser.add_argument("--end-year", type=int, default=None, help="Last observation year (inclusive).")
    parser.add_argument(
        "--refresh-manifest",
        action="store_true",
        help="Refresh complete upstream inventories instead of reusing saved CSV files.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    root = args.output_dir.expanduser().resolve()
    try:
        plan = create_download_plan(
            root,
            argo_source=args.argo_source,
            dataset_version=args.dataset_version,
            start_year=args.start_year,
            end_year=args.end_year,
            refresh=args.refresh_manifest,
        )
        _print_plan(root, plan)
        if args.command == "download":
            download_plan(root, plan)
            print(f"Download complete under {root / 'source'}")
    except (RuntimeError, ValueError) as exc:
        print(f"error: {exc}")
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
