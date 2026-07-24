"""
Render a virtual gel from selected raw read lengths.

The gel uses the same basic signal as the report read-length plot: read length
on the y-axis and DNA mass approximated by read-length weighting. It avoids a
SciPy dependency by histogramming lengths and applying a small Gaussian smooth.
"""

from __future__ import annotations

import math
import os
import tempfile
from dataclasses import dataclass
from collections.abc import Iterable, Mapping
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "altabiotech_mplcache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

DEFAULT_MIN_DISPLAY_BP = 1000
DEFAULT_Y_BINS = 900
DEFAULT_WEIGHT_BY_MASS = True
DEFAULT_MULTIMER_TOLERANCE_FRACTION = 0.15
DEFAULT_MAX_MULTIPLE = 4
DEFAULT_MIN_BAND_READ_COUNT = 5
DEFAULT_MIN_BAND_MASS_FRACTION = 0.02


@dataclass(frozen=True)
class VirtualGelBand:
    center_bp: float
    mass: float
    read_count: int
    multiple: int | None = None


def filtered_display_lengths(lengths: Iterable[int | float], min_display_bp: int) -> np.ndarray:
    values = np.asarray(list(lengths), dtype=float)
    if values.size == 0:
        return values
    values = values[np.isfinite(values)]
    return values[values > min_display_bp]


def _gaussian_kernel(sigma_bins: float) -> np.ndarray:
    sigma_bins = max(float(sigma_bins), 0.5)
    radius = max(1, int(math.ceil(sigma_bins * 4.0)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
    return kernel / kernel.sum()


def resolve_virtual_gel_y_max(length_arrays: list[np.ndarray], y_max: int | None, min_display_bp: int) -> int:
    if y_max is not None:
        if y_max <= 0:
            raise ValueError("y_max must be positive")
        return int(y_max)

    max_values = [float(values.max()) for values in length_arrays if values.size]
    if not max_values:
        return max(2000, min_display_bp * 2)
    rounded = int(math.ceil(max(max_values) / 500.0) * 500)
    return max(rounded, max(2000, min_display_bp * 2))


def _round_axis_top(value: float, min_display_bp: int) -> int:
    rounded = int(math.ceil(value / 500.0) * 500)
    return max(rounded, max(2000, min_display_bp * 2))


def _band_width_bp(center_bp: float) -> float:
    return max(45.0, min(140.0, center_bp * 0.006))


def _weighted_center(lengths: np.ndarray, fallback: float) -> float:
    if lengths.size == 0:
        return fallback
    total = float(lengths.sum())
    if total <= 0:
        return fallback
    return float(np.average(lengths, weights=lengths))


def select_virtual_gel_bands(
    lengths: Iterable[int | float],
    *,
    expected_length_bp: int | float | None = None,
    min_display_bp: int = DEFAULT_MIN_DISPLAY_BP,
    weight_by_mass: bool = DEFAULT_WEIGHT_BY_MASS,
    tolerance_fraction: float = DEFAULT_MULTIMER_TOLERANCE_FRACTION,
    max_multiple: int = DEFAULT_MAX_MULTIPLE,
    min_band_read_count: int = DEFAULT_MIN_BAND_READ_COUNT,
    min_band_mass_fraction: float = DEFAULT_MIN_BAND_MASS_FRACTION,
) -> list[VirtualGelBand]:
    display_lengths = filtered_display_lengths(lengths, min_display_bp)
    if display_lengths.size == 0:
        return []

    weights = display_lengths if weight_by_mass else np.ones_like(display_lengths)
    total_mass = float(weights.sum())
    if total_mass <= 0:
        return []

    if expected_length_bp is not None and expected_length_bp > 0:
        bands = []
        expected = float(expected_length_bp)
        for multiple in range(1, max_multiple + 1):
            expected_center = expected * multiple
            lower = expected_center * (1.0 - tolerance_fraction)
            upper = expected_center * (1.0 + tolerance_fraction)
            in_window = display_lengths[(display_lengths >= lower) & (display_lengths <= upper)]
            if in_window.size < min_band_read_count:
                continue
            band_mass = float(in_window.sum()) if weight_by_mass else float(in_window.size)
            if band_mass / total_mass < min_band_mass_fraction:
                continue
            bands.append(
                VirtualGelBand(
                    center_bp=_weighted_center(in_window, expected_center),
                    mass=band_mass,
                    read_count=int(in_window.size),
                    multiple=multiple,
                )
            )
        return bands

    return _select_density_peaks(
        display_lengths,
        weight_by_mass=weight_by_mass,
        min_band_read_count=min_band_read_count,
        min_band_mass_fraction=min_band_mass_fraction,
    )


def _select_density_peaks(
    lengths: np.ndarray,
    weight_by_mass: bool,
    min_band_read_count: int,
    min_band_mass_fraction: float,
) -> list[VirtualGelBand]:
    if lengths.size == 0:
        return []
    y_max = _round_axis_top(float(lengths.max()), DEFAULT_MIN_DISPLAY_BP)
    bin_width = max(75.0, y_max / 220.0)
    edges = np.arange(0, y_max + bin_width, bin_width)
    weights = lengths if weight_by_mass else None
    density, _ = np.histogram(lengths, bins=edges, weights=weights)
    kernel = _gaussian_kernel(1.5)
    density = np.convolve(density.astype(float), kernel, mode="same")

    total_mass = float((lengths if weight_by_mass else np.ones_like(lengths)).sum())
    if density.size < 3 or total_mass <= 0:
        return []

    bands = []
    peak_threshold = max(float(density.max()) * 0.08, total_mass * 0.01)
    for index in range(1, len(density) - 1):
        if density[index] < density[index - 1] or density[index] < density[index + 1]:
            continue
        if density[index] < peak_threshold:
            continue
        lower = edges[max(0, index - 1)]
        upper = edges[min(len(edges) - 1, index + 2)]
        in_window = lengths[(lengths >= lower) & (lengths < upper)]
        if in_window.size < min_band_read_count:
            continue
        band_mass = float(in_window.sum()) if weight_by_mass else float(in_window.size)
        if band_mass / total_mass < min_band_mass_fraction:
            continue
        bands.append(
            VirtualGelBand(
                center_bp=_weighted_center(in_window, (edges[index] + edges[index + 1]) / 2.0),
                mass=band_mass,
                read_count=int(in_window.size),
            )
        )
    return bands


def resolve_virtual_gel_band_y_max(
    band_groups: list[list[VirtualGelBand]],
    y_max: int | None,
    min_display_bp: int,
) -> int:
    if y_max is not None:
        if y_max <= 0:
            raise ValueError("y_max must be positive")
        return int(y_max)
    max_band = max(
        (band.center_bp + _band_width_bp(band.center_bp) * 6.0 for bands in band_groups for band in bands),
        default=0.0,
    )
    if max_band <= 0:
        return max(2000, min_display_bp * 2)
    return _round_axis_top(max_band * 1.04, min_display_bp)


def _lane_density_from_bands(
    bands: list[VirtualGelBand],
    y_grid: np.ndarray,
) -> np.ndarray:
    density = np.zeros_like(y_grid, dtype=float)
    if not bands:
        return density
    max_mass = max(band.mass for band in bands)
    if max_mass <= 0:
        return density

    for band in bands:
        amplitude = (band.mass / max_mass) ** 0.55
        sigma = _band_width_bp(band.center_bp)
        density += amplitude * np.exp(-0.5 * ((y_grid - band.center_bp) / sigma) ** 2)
    density = np.clip(density, 0.0, 1.0)
    return density


def make_virtual_gel(
    sample_lengths: Mapping[str, Iterable[int | float]],
    out_path: Path,
    title: str,
    *,
    min_display_bp: int = DEFAULT_MIN_DISPLAY_BP,
    y_max: int | None = None,
    y_bins: int = DEFAULT_Y_BINS,
    weight_by_mass: bool = DEFAULT_WEIGHT_BY_MASS,
    expected_lengths: Mapping[str, int | float | None] | None = None,
) -> Path:
    if not sample_lengths:
        raise ValueError("at least one sample is required to render a virtual gel")

    items = [(str(label), lengths) for label, lengths in sample_lengths.items()]
    labels = [label for label, _lengths in items]
    band_groups = [
        select_virtual_gel_bands(
            lengths,
            expected_length_bp=(expected_lengths or {}).get(label),
            min_display_bp=min_display_bp,
            weight_by_mass=weight_by_mass,
        )
        for label, lengths in items
    ]
    if any(band_groups):
        y_max = resolve_virtual_gel_band_y_max(band_groups, y_max, min_display_bp)
    else:
        length_arrays = [filtered_display_lengths(lengths, min_display_bp) for _label, lengths in items]
        y_max = resolve_virtual_gel_y_max(length_arrays, y_max, min_display_bp)

    lane_count = len(labels)
    max_label_len = max(len(label) for label in labels)
    fig_width = max(6.0, min(24.0, 1.15 * lane_count + 2.2))
    fig_height = 7.4 if max_label_len <= 28 else 8.6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    band_cmap = LinearSegmentedColormap.from_list("virtual_gel", ["white", "#202020"])
    lane_width = 0.74
    y_grid = np.linspace(0, y_max, y_bins)
    for index, bands in enumerate(band_groups):
        density = _lane_density_from_bands(bands, y_grid)
        image = np.tile(density[:, None], (1, 48))
        ax.imshow(
            image,
            extent=(index - lane_width / 2.0, index + lane_width / 2.0, 0, y_max),
            origin="lower",
            aspect="auto",
            cmap=band_cmap,
            vmin=0,
            vmax=1,
            interpolation="bicubic",
        )

    ax.set_xlim(-0.55, lane_count - 0.45)
    ax.set_ylim(0, y_max)
    ax.set_xticks(range(lane_count))
    ax.set_xticklabels(labels, rotation=55, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", labelsize=9 if lane_count <= 12 else 8)
    ax.set_ylabel("bp")
    ax.set_title(title, pad=10)
    ax.grid(False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    bottom = 0.28 if max_label_len <= 28 else 0.38
    fig.subplots_adjust(left=0.13, right=0.99, top=0.92, bottom=bottom)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path
