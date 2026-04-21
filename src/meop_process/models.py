from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_VERSION = "MEOP-CTD_yyyy-mm-dd"


@dataclass(frozen=True)
class DiagnosticsDefaults:
    qf: str = "lr1"
    adjusted: bool = True
    parts: tuple[str, ...] = ("tag", "deployment", "overview")


@dataclass(frozen=True)
class BatchDefaults:
    jobs: int = 1
    verbose: bool = False


@dataclass(frozen=True)
class EmailTransportSettings:
    transport: str = "smtp"
    host: str = ""
    port: int = 587
    starttls: bool = True
    use_ssl: bool = False
    username_env: str = ""
    password_env: str = ""
    from_address: str = ""
    sendmail_path: str = "/usr/sbin/sendmail"


@dataclass(frozen=True)
class EmailNotificationSettings:
    enabled: bool = False
    when: str = "always"
    to: tuple[str, ...] = ()
    attach: tuple[str, ...] = ("summary_md",)
    subject_prefix: str = "[MEOP]"
    transport: EmailTransportSettings = field(default_factory=EmailTransportSettings)


@dataclass(frozen=True)
class MeopConfig:
    """Resolved runtime configuration for a pure-Python MEOP checkout or deployment."""

    processdir: Path
    datadir: Path
    publicdir: Path
    pdflatex: str = "pdflatex"
    version: str = DEFAULT_VERSION
    machine: str | None = None
    config_path: Path | None = None
    diagnostics_defaults: DiagnosticsDefaults = field(default_factory=DiagnosticsDefaults)
    batch_defaults: BatchDefaults = field(default_factory=BatchDefaults)
    email_notifications: EmailNotificationSettings = field(default_factory=EmailNotificationSettings)

    @property
    def publicdir_ctd(self) -> Path:
        return self.publicdir / self.version

    @property
    def tablesdir(self) -> Path:
        return self.datadir / "tables"

    @property
    def catalogdir(self) -> Path:
        return self.datadir / "catalog"

    @property
    def data_raw_dir(self) -> Path:
        return self.datadir / "data_raw"

    @property
    def config_files_dir(self) -> Path:
        return self.data_raw_dir / "config_files"

    @property
    def raw_odv_dir(self) -> Path:
        return self.data_raw_dir / "raw_smru_data_odv"

    @property
    def raw_hr_dir(self) -> Path:
        return self.data_raw_dir / "raw_smru_hr_data"

    @property
    def crawl_locdir(self) -> Path:
        return self.data_raw_dir / "crawl_locations"

    @property
    def cls_locdir(self) -> Path:
        return self.data_raw_dir / "smooth_cls_locations"

    @property
    def final_dataset_dir(self) -> Path:
        return self.datadir / "data_prof"

    @property
    def trajectory_dataset_dir(self) -> Path:
        return self.datadir / "data_traj"

    @property
    def mapsdir(self) -> Path:
        return self.datadir / "maps"

    @property
    def plotdir(self) -> Path:
        return self.datadir / "plots_by_tags"

    @property
    def plots_by_deployment_dir(self) -> Path:
        return self.datadir / "plots_by_deployments"

    @property
    def plots_overview_dir(self) -> Path:
        return self.datadir / "plots_overview"


@dataclass(frozen=True)
class Selection:
    """Target deployment and/or tag selection."""

    deployment: str = ""
    smru_name: str = ""

    def normalized(self) -> "Selection":
        if self.smru_name and not self.deployment:
            return Selection(
                deployment=self.smru_name.split("-")[0],
                smru_name=self.smru_name,
            )
        return self

    @property
    def label(self) -> str:
        if self.smru_name:
            return self.smru_name
        return self.deployment


@dataclass(frozen=True)
class DeploymentRecord:
    """One deployment entry resolved from list_deployment.csv and mirrored JSON metadata."""

    deployment_code: str
    pi_code: str = ""
    country: str = ""
    process: str = ""
    public: str = ""
    description: str = ""
    gts: str = ""
    start_date: str = ""
    end_date: str = ""
    task_done: str = ""
    source: str = "catalog"
    row_name: str = ""
    extra: dict[str, str] = field(default_factory=dict)

    @property
    def pi(self) -> str:
        return self.pi_code.strip()

    @property
    def nation(self) -> str:
        return self.country

    def as_dict(self) -> dict[str, Any]:
        return {
            "deployment_code": self.deployment_code,
            "pi_code": self.pi_code,
            "country": self.country,
            "process": self.process,
            "public": self.public,
            "description": self.description,
            "gts": self.gts,
            "start_date": self.start_date,
            "end_date": self.end_date,
            "task_done": self.task_done,
            "source": self.source,
            "row_name": self.row_name or self.deployment_code,
            **self.extra,
        }


@dataclass(frozen=True)
class DeploymentInfo:
    """Python-native replacement for the legacy ``info_deployment`` structure."""

    selection: Selection
    record: DeploymentRecord | None
    invalid_code: bool
    directory: Path
    raw_input_dir: Path
    raw_input_zip: Path
    raw_working_text: Path
    raw_working_ctd_text: Path | None = None
    raw_working_fl_text: Path | None = None
    raw_working_fcell: Path | None = None
    raw_smru_names: tuple[str, ...] = ()
    raw_profile_count_by_smru: dict[str, int] = field(default_factory=dict)
    list_smru_name: tuple[str, ...] = ()
    list_tag_lr0: tuple[Path, ...] = ()
    list_tag_lr1: tuple[Path, ...] = ()
    list_tag_hr0: tuple[Path, ...] = ()
    list_tag_hr1: tuple[Path, ...] = ()
    list_tag_hr2: tuple[Path, ...] = ()
    list_tag_fr0: tuple[Path, ...] = ()
    list_tag_fr1: tuple[Path, ...] = ()
    known_platform_codes: tuple[str, ...] = ()
    hr_platform_codes: tuple[str, ...] = ()

    @property
    def EXP(self) -> str:
        return self.selection.deployment

    @property
    def PI(self) -> str:
        return self.record.pi if self.record is not None else ""

    @property
    def NATION(self) -> str:
        return self.record.nation if self.record is not None else ""

    @property
    def nomfic(self) -> str:
        if self.selection.deployment in {"ct3", "ct7", "ct11", "wd3"}:
            return f"{self.selection.deployment}_fcell.mat"
        return f"{self.selection.deployment}_ODV.txt"

    @property
    def process(self) -> str:
        return self.record.process if self.record is not None else ""

    @property
    def public(self) -> str:
        return self.record.public if self.record is not None else ""

    @property
    def requested_smru_name(self) -> str:
        return self.selection.smru_name

    @property
    def deployment_exists(self) -> bool:
        return not self.invalid_code and self.record is not None

    def as_dict(self) -> dict[str, Any]:
        return {
            "deployment": self.selection.deployment,
            "smru_name": self.selection.smru_name,
            "invalid_code": self.invalid_code,
            "pi": self.PI,
            "nation": self.NATION,
            "nomfic": self.nomfic,
            "directory": str(self.directory),
            "raw_input_dir": str(self.raw_input_dir),
            "raw_input_zip": str(self.raw_input_zip),
            "raw_working_text": str(self.raw_working_text),
            "raw_working_ctd_text": str(self.raw_working_ctd_text) if self.raw_working_ctd_text else "",
            "raw_working_fl_text": str(self.raw_working_fl_text) if self.raw_working_fl_text else "",
            "raw_working_fcell": str(self.raw_working_fcell) if self.raw_working_fcell else "",
            "raw_smru_names": list(self.raw_smru_names),
            "raw_profile_count_by_smru": dict(self.raw_profile_count_by_smru),
            "list_smru_name": list(self.list_smru_name),
            "list_tag_lr0": [str(path) for path in self.list_tag_lr0],
            "list_tag_lr1": [str(path) for path in self.list_tag_lr1],
            "list_tag_hr0": [str(path) for path in self.list_tag_hr0],
            "list_tag_hr1": [str(path) for path in self.list_tag_hr1],
            "list_tag_hr2": [str(path) for path in self.list_tag_hr2],
            "list_tag_fr0": [str(path) for path in self.list_tag_fr0],
            "list_tag_fr1": [str(path) for path in self.list_tag_fr1],
            "known_platform_codes": list(self.known_platform_codes),
            "hr_platform_codes": list(self.hr_platform_codes),
        }
