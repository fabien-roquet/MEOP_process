from __future__ import annotations

from meop_process.config.sync import sync_external_config_files


def test_sync_external_config_files_updates_and_backs_up_previous_copy(meop_config) -> None:
    source_dir = meop_config.datadir / "external_config_files"
    source_dir.mkdir(parents=True, exist_ok=True)
    source = source_dir / "deployment3.json"
    source.write_text('{"version": 2}', encoding="utf-8")

    destination_dir = meop_config.config_files_dir
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / "deployment3.json"
    destination.write_text('{"version": 1}', encoding="utf-8")

    updated = sync_external_config_files(meop_config, timestamp="20260318")

    assert destination in updated
    assert destination.read_text(encoding="utf-8") == '{"version": 2}'
    assert (destination_dir / "deployment3_20260318.json").read_text(encoding="utf-8") == '{"version": 1}'
