#!/usr/bin/env python3
"""Plot the latitude–P/Ps statistics written by lonlat_to_3D.py.

Requires NumPy, Matplotlib, and VTK. Edit the settings below or run with --help.
"""
import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLImageDataReader


# Editable defaults. INPUT_DIR contains the *_statistics_zonal_*.vti files.
RUN = "SimpleJ7Z30"
SEASON = "summer"
INPUT_DIR = Path(".")
OUTPUT_DIR = Path(".")
DPI = 160
SHOW_FIGURES = False
LATITUDE_LIMITS = (-90, 90)
PRESSURE_LIMITS = (1, 0)  # Surface at the bottom, top of atmosphere at the top.
MEAN_CMAP = "jet"
SPREAD_CMAP = "viridis"
N_CONTOUR_BANDS = 32

# Panel order, colorbar label, mean color limits, and colorbar tick step.
# These match the ParaView reference. Summer has a different variance scale.
FIELDS = (
    ("Temperature", "Temperature [K]", (140, 340), 20),
    ("Velocity_Zonal", r"Zonal velocity [m s$^{-1}$]", (-40, 40), 5),
    ("EddyKineticEnergy", r"Eddy Kinetic Energy [m$^2$ s$^{-2}$]", (-600, 600), 100),
    ("EddyHeatFlux", r"Eddy Heat Flux [K m s$^{-1}$]", (-40, 40), 5),
    ("TemperatureVariance", r"Temperature Variance [K$^2$]", (-80, 80), 10),
    ("EddyMomentumFlux", r"Eddy Momentum Flux [m$^2$ s$^{-2}$]", (-100, 100), 20),
)
MEAN_SCALE_OVERRIDES = {"summer": {"TemperatureVariance": ((-80, 80), 10)}}

# Layout is one six-panel figure per statistic. For fixed spread scales across
# runs, add a maximum keyed by (suffix, field), e.g. ("_TimeStdDev", "Temperature"): 10.
SPREAD_MAX_OVERRIDES = {}
STATISTICS = (
    ("", "mean", "Climatological mean"),
    ("_TimeStdDev", "stddev", "Temporal standard deviation"),
    ("_CI95HalfWidth", "ci95_half_width", "95% confidence half-width of the mean"),
)


def input_path(explicit_path, directory, run, season):
    if explicit_path is not None:
        if not explicit_path.is_file():
            raise FileNotFoundError(f"Input VTI file not found: {explicit_path}")
        return explicit_path
    path = directory / f"{run}_statistics_zonal_{season}.vti"
    if path.is_file():
        return path
    # Prefer summer by default, with a fallback for runs that only have all data.
    all_path = directory / f"{run}_statistics_zonal_all.vti"
    if season == "summer" and all_path.is_file():
        print(f"{path.name} is unavailable. Using {all_path.name}.")
        return all_path
    raise FileNotFoundError(f"Input VTI file not found: {path}")


def read_zonal_vti(path):
    reader = vtkXMLImageDataReader()
    if not reader.CanReadFile(str(path)):
        raise ValueError(f"Cannot read VTI file: {path}")
    reader.SetFileName(str(path))
    reader.Update()
    image = reader.GetOutput()
    nx, nlat, nz = image.GetDimensions()
    if nx != 1 or nlat < 2 or nz < 2:
        raise ValueError(f"{path.name} must have zonal VTI dimensions (1, nlat, nz).")

    origin = image.GetOrigin()
    spacing = image.GetSpacing()
    extent = image.GetExtent()
    lat = origin[1] + spacing[1] * (extent[2] + np.arange(nlat))
    pressure = origin[2] + spacing[2] * (extent[4] + np.arange(nz))

    point_data = image.GetPointData()
    names = [field + suffix for field, _, _, _ in FIELDS for suffix, _, _ in STATISTICS]
    arrays = {}
    for name in names:
        vtk_array = point_data.GetArray(name)
        if vtk_array is None:
            continue
        if vtk_array.GetNumberOfComponents() != 1:
            raise ValueError(f"Expected a scalar field for {name!r}.")
        # VTK stores x fastest. With nx=1, latitude varies fastest.
        arrays[name] = vtk_to_numpy(vtk_array).reshape(nz, nlat).copy()
    if not arrays:
        raise ValueError(f"{path.name} contains none of the expected zonal statistics arrays.")
    return lat, pressure, arrays


def rounded_maximum(value):
    """Round up to a readable scale, keeping sub-unit fields on a decimal scale."""
    if value <= 0:
        return 1
    step = 10.0 ** (math.floor(math.log10(value)) - 1)
    if value >= 1:
        step = max(1, step)
    return math.ceil(value / step) * step


def write_peak_confidence(path, arrays):
    """Write one row per field at the largest finite absolute mean."""
    with path.open("w", encoding="utf-8") as report:
        print("Pointwise 95% confidence at the maximum absolute climatological mean", file=report)
        print("All values are evaluated at the grid point where |M| is largest.", file=report)
        print(f"{'Field':<21} {'Mean M':>12} {'Half-width H':>12} "
              f"{'100H/|M|':>10} {'100StdDev/|M|':>15}", file=report)
        for field, _, _, _ in FIELDS:
            mean = arrays.get(field)
            if mean is None or not np.isfinite(mean).any():
                print(f"{field:<21} no finite mean values", file=report)
                continue

            magnitude = np.where(np.isfinite(mean), np.abs(mean), -np.inf)
            k, j = np.unravel_index(np.argmax(magnitude), mean.shape)
            value = mean[k, j]
            half_widths = arrays.get(field + "_CI95HalfWidth")
            half_width = half_widths[k, j] if half_widths is not None else np.nan
            std_devs = arrays.get(field + "_TimeStdDev")
            std_dev = std_devs[k, j] if std_devs is not None else np.nan
            width_text = f"{half_width:.5g}" if np.isfinite(half_width) else "N/A"
            relative = f"{100 * half_width / abs(value):.2f}%" if (
                np.isfinite(half_width) and value != 0) else "N/A"
            relative_std_dev = f"{100 * std_dev / abs(value):.2f}%" if (
                np.isfinite(std_dev) and value != 0) else "N/A"
            print(f"{field:<21} {value:12.5g} {width_text:>12} "
                  f"{relative:>10} {relative_std_dev:>15}", file=report)


def make_figure(lat, pressure, arrays, season, suffix):
    fig, axes = plt.subplots(2, 3, figsize=(17, 7.7), layout="constrained")

    for ax, (field, label, limits, tick_step) in zip(axes.flat, FIELDS):
        name = field + suffix
        values = arrays.get(name)
        ax.set(xlabel=r"Latitude [$^\circ$]", ylabel=r"$P/P_s$")
        ax.yaxis.label.set_rotation(0)
        ax.yaxis.label.set_horizontalalignment("right")
        ax.set_xticks(np.arange(-90, 91, 30))
        ax.set_yticks(np.linspace(0, 1, 5), ["0", "0.25", "0.5", "0.75", "1"])
        ax.set(xlim=LATITUDE_LIMITS, ylim=PRESSURE_LIMITS)
        finite = values[np.isfinite(values)] if values is not None else np.array([])
        if not finite.size:
            message = "Array missing from VTI" if values is None else "No finite values"
            ax.text(0.5, 0.5, f"{field}\n{message}", ha="center", va="center", transform=ax.transAxes)
            continue

        if suffix:
            maximum = SPREAD_MAX_OVERRIDES.get((suffix, field), rounded_maximum(float(finite.max())))
            lo, hi = 0, maximum
            ticks = np.linspace(lo, hi, 6)
            cmap = SPREAD_CMAP
        else:
            (lo, hi), tick_step = MEAN_SCALE_OVERRIDES.get(season, {}).get(
                field, (limits, tick_step))
            ticks = np.arange(lo, hi + 0.5 * tick_step, tick_step)
            cmap = MEAN_CMAP

        # Saturate values outside the display range, as in the ParaView panels.
        levels = np.linspace(lo, hi, N_CONTOUR_BANDS + 1)
        display_values = np.clip(np.ma.masked_invalid(values), lo, hi)
        contours = ax.contourf(lat, pressure, display_values, levels=levels, cmap=cmap)
        colorbar = fig.colorbar(contours, ax=ax, ticks=ticks, pad=0.025, fraction=0.05,
                               format=FormatStrFormatter("%g"))
        colorbar.set_label(label)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vti", type=Path, nargs="?", help="explicit input VTI path")
    parser.add_argument("--run", default=RUN, help="run prefix")
    parser.add_argument("--season", default=SEASON, choices=("spring", "summer", "fall", "winter", "all"))
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--show", action="store_true", default=SHOW_FIGURES)
    args = parser.parse_args()

    try:
        path = input_path(args.vti, args.input_dir, args.run, args.season)
        lat, pressure, arrays = read_zonal_vti(path)
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))

    prefix, separator, season = path.stem.rpartition("_statistics_zonal_")
    if not separator:
        prefix, season = path.stem, args.season
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / f"{prefix}_climatology_{season}_peak_confidence.txt"
    write_peak_confidence(report_path, arrays)
    print(report_path)
    for suffix, tag, title in STATISTICS:
        expected = [field + suffix for field, _, _, _ in FIELDS]
        missing = [name for name in expected if name not in arrays]
        if len(missing) == len(expected):
            print(f"Skipping {title.lower()}: no matching arrays in {path.name}.", file=sys.stderr)
            continue
        if missing:
            print(f"{path.name} is missing: {', '.join(missing)}.", file=sys.stderr)
        fig = make_figure(lat, pressure, arrays, season, suffix)
        output = args.output_dir / f"{prefix}_climatology_{season}_{tag}.png"
        fig.savefig(output, dpi=DPI)
        print(output)
        if not args.show:
            plt.close(fig)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
