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


def _lane_density(
    lengths: np.ndarray,
    y_max: int,
    y_bins: int,
    weight_by_mass: bool,
) -> np.ndarray:
    if lengths.size == 0:
        return np.zeros(y_bins, dtype=float)

    edges = np.linspace(0, y_max, y_bins + 1)
    weights = lengths if weight_by_mass else None
    density, _ = np.histogram(lengths, bins=edges, weights=weights)

    bin_width = y_max / y_bins
    sigma_bp = max(50.0, y_max / 160.0)
    kernel = _gaussian_kernel(sigma_bp / bin_width)
    density = np.convolve(density.astype(float), kernel, mode="same")

    max_density = density.max()
    if max_density > 0:
        density = density / max_density
        density = np.power(density, 0.6)
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
) -> Path:
    if not sample_lengths:
        raise ValueError("at least one sample is required to render a virtual gel")

    items = [(str(label), lengths) for label, lengths in sample_lengths.items()]
    labels = [label for label, _lengths in items]
    length_arrays = [filtered_display_lengths(lengths, min_display_bp) for _label, lengths in items]
    y_max = resolve_virtual_gel_y_max(length_arrays, y_max, min_display_bp)

    lane_count = len(labels)
    max_label_len = max(len(label) for label in labels)
    fig_width = max(6.0, min(24.0, 1.15 * lane_count + 2.2))
    fig_height = 7.4 if max_label_len <= 28 else 8.6
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    band_cmap = LinearSegmentedColormap.from_list("virtual_gel", ["white", "#202020"])
    lane_width = 0.74
    for index, lengths in enumerate(length_arrays):
        density = _lane_density(
            lengths=lengths,
            y_max=y_max,
            y_bins=y_bins,
            weight_by_mass=weight_by_mass,
        )
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
    fig.subplots_adjust(left=0.09, right=0.99, top=0.92, bottom=bottom)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path
