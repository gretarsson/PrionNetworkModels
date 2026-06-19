#!/usr/bin/env python3
"""Utilities for exporting manuscript panels into per-figure folders.

Figures 2-5 were manually composed before this repository had a clean paper
reproduction layer. These helpers make the current composed scientific panels
available as independent files while keeping explicit provenance notes for the
inference bundles that still need to be translated into `runs/`.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANUSCRIPT_FIGURES = Path(
    "/Users/gretarsson/Desktop/nonlinear_synuclein_manuscript/figures"
)


@dataclass(frozen=True)
class Crop:
    name: str
    source: str
    x: int
    y: int
    width: int
    height: int
    description: str
    note: str = ""


def manuscript_figure_dir() -> Path:
    return Path(
        __import__("os").environ.get(
            "MANUSCRIPT_FIGURES_DIR", str(DEFAULT_MANUSCRIPT_FIGURES)
        )
    )


def crop_with_sips(source: Path, crop: Crop, out_png: Path) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as img:
        box = (crop.x, crop.y, crop.x + crop.width, crop.y + crop.height)
        cropped = img.crop(box)
        cropped.save(out_png)


def image_size(path: Path) -> tuple[int, int]:
    proc = subprocess.run(
        ["sips", "-g", "pixelWidth", "-g", "pixelHeight", str(path)],
        check=True,
        text=True,
        capture_output=True,
    )
    width = height = None
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith("pixelWidth:"):
            width = int(line.split(":", 1)[1])
        elif line.startswith("pixelHeight:"):
            height = int(line.split(":", 1)[1])
    if width is None or height is None:
        raise RuntimeError(f"Could not read image size for {path}")
    return width, height


def png_to_pdf(png: Path, pdf: Path) -> None:
    with Image.open(png) as img:
        if img.mode == "RGBA":
            img = img.convert("RGB")
        img.save(pdf)


def export_crops(figure: str, crops: list[Crop], *, dependencies: list[dict]) -> Path:
    src_root = manuscript_figure_dir()
    out_dir = ROOT / "paper-rf" / "figures" / figure
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "figure": figure,
        "source_directory": str(src_root),
        "panels": [],
        "run_bundle_dependencies": dependencies,
    }

    for crop in crops:
        source = src_root / crop.source
        if not source.exists():
            raise FileNotFoundError(f"Missing source figure: {source}")
        png = out_dir / f"{crop.name}.png"
        pdf = out_dir / f"{crop.name}.pdf"
        crop_with_sips(source, crop, png)
        png_to_pdf(png, pdf)
        manifest["panels"].append(
            {
                "name": crop.name,
                "source": str(source),
                "description": crop.description,
                "note": crop.note,
                "png": str(png),
                "pdf": str(pdf),
            }
        )

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return out_dir


def write_missing_requirements(path: Path, requirements: list[str]) -> None:
    text = ["# Missing Run Bundles", ""]
    if requirements:
        text.extend(f"- {item}" for item in requirements)
    else:
        text.append("None.")
    path.write_text("\n".join(text) + "\n")
