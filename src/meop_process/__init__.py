"""Pure-Python package entrypoint for the MEOP processing workflow."""

from .api import (
    apply_adjustments,
    bootstrap_data_store,
    compare_reference_outputs,
    create_fr0,
    create_hr2,
    describe_runtime_data_layout,
    generate_diagnostics,
    import_raw_data,
    load_info_deployment,
    process_tags,
    run_all_deployments,
    update_config_files,
    update_metadata,
    update_metadata_summaries,
    validate_runtime_data_layout,
)

__all__ = [
    "apply_adjustments",
    "bootstrap_data_store",
    "compare_reference_outputs",
    "create_fr0",
    "create_hr2",
    "describe_runtime_data_layout",
    "generate_diagnostics",
    "import_raw_data",
    "load_info_deployment",
    "process_tags",
    "run_all_deployments",
    "update_config_files",
    "update_metadata",
    "update_metadata_summaries",
    "validate_runtime_data_layout",
]

__version__ = "0.1.0"
