"""Generate a minimal static HTML site from diagnostics plot directories."""
from __future__ import annotations

import os
import textwrap
from pathlib import Path

# ---------------------------------------------------------------------------
# HTML skeleton
# ---------------------------------------------------------------------------

_CSS = textwrap.dedent("""\
    body { font-family: sans-serif; margin: 0; padding: 1rem 2rem; background: #f8f9fa; color: #333; }
    h1, h2 { color: #2c5f8a; }
    nav { margin-bottom: 1.5rem; font-size: 0.95rem; }
    nav a { margin-right: 1rem; color: #2c5f8a; text-decoration: none; }
    nav a:hover { text-decoration: underline; }
    .plots { display: flex; flex-wrap: wrap; gap: 0.75rem; margin-bottom: 1.5rem; }
    .plots img { max-width: 480px; height: auto; border: 1px solid #ccc; background: #fff; }
    ul.tag-list { columns: 3; list-style: none; padding: 0; }
    ul.tag-list li { margin-bottom: 0.4rem; }
    ul.tag-list a { color: #2c5f8a; text-decoration: none; }
    ul.tag-list a:hover { text-decoration: underline; }
    ul.dep-list { columns: 4; list-style: none; padding: 0; }
    ul.dep-list li { margin-bottom: 0.4rem; }
    ul.dep-list a { color: #2c5f8a; text-decoration: none; }
    ul.dep-list a:hover { text-decoration: underline; }
    pre.info { background: #f0f0f0; padding: 1rem; white-space: pre-wrap; font-size: 0.9rem; }
""")

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
{css}</style>
</head>
<body>
{nav}<h1>{heading}</h1>
{body}
</body>
</html>
"""


def _write_html(path: Path, title: str, heading: str, body: str, nav_links: list[tuple[str, Path]] = ()) -> Path:
    """Render an HTML page and write it to *path*."""
    page_dir = path.parent
    nav_html = ""
    if nav_links:
        anchors = " ".join(
            f'<a href="{_relpath(page_dir, target)}">{label}</a>'
            for label, target in nav_links
        )
        nav_html = f"<nav>{anchors}</nav>\n"

    html = _HTML_TEMPLATE.format(title=title, heading=heading, css=_CSS, nav=nav_html, body=body)
    path.write_text(html, encoding="utf-8")
    return path


def _relpath(from_dir: Path, to: Path) -> str:
    """Relative URL from *from_dir* to *to* using forward slashes."""
    return os.path.relpath(to, from_dir).replace(os.sep, "/")


def _images_html(images: list[Path], page_dir: Path) -> str:
    """Return a <div class="plots"> block with lazy-loaded <img> tags."""
    if not images:
        return ""
    items = "\n".join(
        f'  <img src="{_relpath(page_dir, img)}" alt="{img.stem}" loading="lazy">'
        for img in sorted(images)
    )
    return f'<div class="plots">\n{items}\n</div>'


# ---------------------------------------------------------------------------
# SMRU-name extraction helper
# ---------------------------------------------------------------------------

def _smru_from_stem(stem: str) -> str | None:
    """Extract the SMRU platform code from a plot file stem.

    Plot filenames follow ``{smru}_{qf}_{type}_{suffix}``.  The QF token is
    ``lr1``, ``lr2``, or ``hr2``.  We split on the first ``_lr`` / ``_hr``
    occurrence.
    """
    for sep in ("_lr", "_hr"):
        idx = stem.find(sep)
        if idx > 0:
            return stem[:idx]
    return None


# ---------------------------------------------------------------------------
# Page builders
# ---------------------------------------------------------------------------

def build_overview_page(
    plots_overview: Path,
    plots_by_deployments: Path,
    *,
    rebuild: bool = False,
) -> Path | None:
    """Build ``plots_overview/index.html`` listing all deployments."""
    if not plots_overview.is_dir():
        return None
    out = plots_overview / "index.html"
    if out.exists() and not rebuild:
        return out

    images = sorted(plots_overview.glob("*.png"))
    plots_html = _images_html(images, out.parent)

    dep_names = sorted(d.name for d in plots_by_deployments.iterdir() if d.is_dir()) \
        if plots_by_deployments.is_dir() else []
    dep_items = "\n".join(
        f'  <li><a href="{_relpath(out.parent, plots_by_deployments / d / "index.html")}">{d}</a></li>'
        for d in dep_names
    )
    dep_section = f'<h2>Deployments ({len(dep_names)})</h2>\n<ul class="dep-list">\n{dep_items}\n</ul>' \
        if dep_names else ""

    body = f"{plots_html}\n{dep_section}"
    return _write_html(out, "MEOP-CTD Overview", "MEOP-CTD Overview", body)


def build_deployment_pages(
    plots_by_deployments: Path,
    plots_by_tags: Path,
    *,
    rebuild: bool = False,
) -> list[Path]:
    """Build ``plots_by_deployments/{dep}/index.html`` for every deployment."""
    if not plots_by_deployments.is_dir():
        return []

    written: list[Path] = []
    overview_index = plots_by_deployments.parent / "plots_overview" / "index.html"

    for dep_dir in sorted(plots_by_deployments.iterdir()):
        if not dep_dir.is_dir():
            continue
        out = dep_dir / "index.html"
        if out.exists() and not rebuild:
            written.append(out)
            continue

        dep = dep_dir.name
        images = [p for p in sorted(dep_dir.glob("*.png"))]
        plots_html = _images_html(images, out.parent)

        # Discover tags for this deployment
        tag_dir = plots_by_tags / dep
        smru_names: list[str] = []
        if tag_dir.is_dir():
            seen: set[str] = set()
            for img in tag_dir.glob("*.png"):
                smru = _smru_from_stem(img.stem)
                if smru:
                    seen.add(smru)
            smru_names = sorted(seen)

        tag_items = "\n".join(
            f'  <li><a href="{_relpath(out.parent, plots_by_tags / dep / (smru + ".html"))}">{smru}</a></li>'
            for smru in smru_names
        )
        tag_section = f'<h2>Tags ({len(smru_names)})</h2>\n<ul class="tag-list">\n{tag_items}\n</ul>' \
            if smru_names else ""

        body = f"{plots_html}\n{tag_section}"
        nav: list[tuple[str, Path]] = []
        if overview_index.exists() or True:  # always include; file may be created later
            nav.append(("Overview", overview_index))

        _write_html(out, f"Deployment {dep}", f"Deployment: {dep}", body, nav_links=nav)
        written.append(out)

    return written


def build_tag_pages(
    plots_by_tags: Path,
    plots_by_deployments: Path,
    *,
    rebuild: bool = False,
) -> list[Path]:
    """Build ``plots_by_tags/{dep}/{smru}.html`` for every tag."""
    if not plots_by_tags.is_dir():
        return []

    written: list[Path] = []
    overview_index = plots_by_tags.parent / "plots_overview" / "index.html"

    for dep_dir in sorted(plots_by_tags.iterdir()):
        if not dep_dir.is_dir():
            continue
        dep = dep_dir.name

        # Group PNG files and info TXT by SMRU name
        smru_images: dict[str, list[Path]] = {}
        smru_info: dict[str, Path] = {}

        for img in dep_dir.glob("*.png"):
            smru = _smru_from_stem(img.stem)
            if smru:
                smru_images.setdefault(smru, []).append(img)

        for txt in dep_dir.glob("*.txt"):
            smru = _smru_from_stem(txt.stem)
            if smru:
                smru_info[smru] = txt

        dep_index = plots_by_deployments / dep / "index.html"

        for smru in sorted(smru_images):
            out = dep_dir / f"{smru}.html"
            if out.exists() and not rebuild:
                written.append(out)
                continue

            images = smru_images[smru]
            plots_html = _images_html(images, out.parent)

            info_html = ""
            if smru in smru_info:
                try:
                    text = smru_info[smru].read_text(encoding="utf-8")
                    info_html = f'<h2>Info</h2>\n<pre class="info">{text}</pre>'
                except OSError:
                    pass

            nav: list[tuple[str, Path]] = []
            nav.append(("Overview", overview_index))
            nav.append((dep, dep_index))

            body = f"{plots_html}\n{info_html}"
            _write_html(out, f"Tag {smru}", f"Tag: {smru}", body, nav_links=nav)
            written.append(out)

    return written


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_site(
    plots_by_tags: Path,
    plots_by_deployments: Path,
    plots_overview: Path,
    *,
    rebuild: bool = False,
    verbose: bool = False,
) -> list[Path]:
    """Generate all static HTML pages for the diagnostics site.

    Creates:
    - ``plots_overview/index.html``
    - ``plots_by_deployments/{dep}/index.html`` for each deployment
    - ``plots_by_tags/{dep}/{smru}.html`` for each tag

    Parameters
    ----------
    plots_by_tags:
        Directory containing per-deployment subdirectories of tag plots.
    plots_by_deployments:
        Directory containing per-deployment subdirectories of deployment plots.
    plots_overview:
        Directory containing cross-deployment overview plots.
    rebuild:
        Overwrite existing HTML files.
    verbose:
        Print progress messages to stdout.

    Returns
    -------
    list[Path]
        Paths to HTML files written (or already present when ``rebuild=False``).
    """
    written: list[Path] = []

    overview = build_overview_page(plots_overview, plots_by_deployments, rebuild=rebuild)
    if overview is not None:
        written.append(overview)
        if verbose:
            print(f"  site overview : {overview}")

    dep_pages = build_deployment_pages(plots_by_deployments, plots_by_tags, rebuild=rebuild)
    written.extend(dep_pages)
    if verbose and dep_pages:
        print(f"  deployment pages : {len(dep_pages)}")

    tag_pages = build_tag_pages(plots_by_tags, plots_by_deployments, rebuild=rebuild)
    written.extend(tag_pages)
    if verbose and tag_pages:
        print(f"  tag pages : {len(tag_pages)}")

    return written
