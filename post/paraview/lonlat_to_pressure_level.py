#!/usr/bin/env python3
# Usage examples:
#   python lonlat_to_pressure_level_mean_native_griddata.py SimpleJ5Z30 5 5 30 30 163
#   python lonlat_to_pressure_level_mean_native_griddata.py SimpleJ5Z30 5 5 30 30 163 --pressure-hpa 500
#   srun -n 30 python -m mpi4py lonlat_to_pressure_level_mean_native_griddata.py SimpleJ5Z30 5 5 30 30 163 --pressure-hpa 500
#
# Generates one time-averaged 2D lon-lat *.vti file with fields interpolated
# directly from native uneven pressure columns to a specified pressure level.
# Default pressure level is 500 hPa.
#
# Algorithm:
#   1. Convert each time to lon-lat VTP layers with run_xyz_layers().
#   2. For each native horizontal cell, compute P_k = (P/Ps)_k * Ps_k.
#   3. Interpolate each cell-data field vertically in the native pressure column.
#   4. Interpolate the resulting native cell-centred 2D field horizontally to a
#      regular lon-lat image using scipy.interpolate.griddata.
#   5. Accumulate a pointwise running mean over all requested times.

import os
import sys
import argparse
import textwrap
from contextlib import suppress

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from scipy.interpolate import griddata

from utilities import *
from xyz2lonlat_layers import run_xyz_layers
from mpi4py import MPI


def read_vtp(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    out = vtk.vtkPolyData()
    out.DeepCopy(reader.GetOutput())
    return out


def pressure_target_in_data_units(ps_values, pressure_hpa):
    """Return target pressure in the same units as Ps."""
    finite_ps = ps_values[np.isfinite(ps_values)]
    if finite_ps.size == 0:
        raise RuntimeError("Cannot determine Ps units: no finite Ps values found.")

    # Typical surface pressure is about 100000 Pa or 1000 hPa.
    # If median Ps is below 2000, assume hPa; otherwise assume Pa.
    if np.nanmedian(finite_ps) < 2000.0:
        return pressure_hpa
    return pressure_hpa * 100.0


def interp_column(pcol, vcol, target_pressure):
    """Piecewise-linear interpolation for one native uneven pressure column."""
    pcol = np.asarray(pcol, dtype=np.float64)
    vcol = np.asarray(vcol, dtype=np.float64)

    if vcol.ndim == 1:
        ok = np.isfinite(pcol) & np.isfinite(vcol)
        if np.count_nonzero(ok) < 2:
            return np.nan
        p = pcol[ok]
        v = vcol[ok]
        order = np.argsort(p)
        p = p[order]
        v = v[order]
        p, unique_idx = np.unique(p, return_index=True)
        v = v[unique_idx]
        if p.size < 2 or target_pressure < p[0] or target_pressure > p[-1]:
            return np.nan
        return np.interp(target_pressure, p, v)

    # Multi-component array, shape (nz, ncomp).
    out = np.full(vcol.shape[1], np.nan, dtype=np.float64)
    for c in range(vcol.shape[1]):
        out[c] = interp_column(pcol, vcol[:, c], target_pressure)
    return out


def cell_centres_lonlat(polydata):
    """Return native cell centres as (lon, lat) from the VTP horizontal mesh."""
    coords = vtk_to_numpy(polydata.GetPoints().GetData())
    polys = vtk_to_numpy(polydata.GetPolys().GetData())
    n_cells = polydata.GetNumberOfCells()

    centres = np.empty((n_cells, 2), dtype=np.float64)
    pos = 0
    for icell in range(n_cells):
        n = int(polys[pos])
        ids = polys[pos + 1:pos + 1 + n]
        xy = coords[ids, :2]
        centres[icell, 0] = np.mean(xy[:, 0])
        centres[icell, 1] = np.mean(xy[:, 1])
        pos += n + 1
    return centres


def make_native_pressure_level_fields(vtp_series, pressure_hpa):
    """
    Return native cell centres and pressure-level fields.

    fields[name] is shape (n_cells,) for scalar arrays or (n_cells,ncomp) for
    vector arrays.
    """
    layers = [read_vtp(f) for f in vtp_series]
    first = layers[0]

    n_cells = first.GetNumberOfCells()
    n_points = first.GetNumberOfPoints()

    for k, ds in enumerate(layers):
        if ds.GetNumberOfCells() != n_cells:
            raise RuntimeError(f"Layer {k + 1} has {ds.GetNumberOfCells()} cells; expected {n_cells}.")
        if ds.GetNumberOfPoints() != n_points:
            raise RuntimeError(f"Layer {k + 1} has {ds.GetNumberOfPoints()} points; expected {n_points}.")

    cell_data0 = first.GetCellData()
    array_names = [cell_data0.GetArrayName(i) for i in range(cell_data0.GetNumberOfArrays())]

    if "P/Ps" not in array_names or "Ps" not in array_names:
        raise RuntimeError("Pressure-level interpolation requires cell-data arrays named 'P/Ps' and 'Ps'.")

    stacks = {}
    for name in array_names:
        arrs = []
        for ds in layers:
            arr = ds.GetCellData().GetArray(name)
            if arr is None:
                raise RuntimeError(f"Array '{name}' is missing from one or more layers.")
            arrs.append(vtk_to_numpy(arr))
        stacks[name] = np.stack(arrs, axis=0)

    p_over_ps = np.asarray(stacks["P/Ps"], dtype=np.float64)
    ps = np.asarray(stacks["Ps"], dtype=np.float64)
    pressure = p_over_ps * ps
    target_pressure = pressure_target_in_data_units(ps, pressure_hpa)

    centres = cell_centres_lonlat(first)
    fields = {}

    for name in array_names:
        data = stacks[name]
        sample = data[0]

        if sample.ndim == 1:
            level = np.full(n_cells, np.nan, dtype=np.float64)
        else:
            level = np.full((n_cells, sample.shape[1]), np.nan, dtype=np.float64)

        for icell in range(n_cells):
            if sample.ndim == 1:
                level[icell] = interp_column(pressure[:, icell], data[:, icell], target_pressure)
            else:
                level[icell, :] = interp_column(pressure[:, icell], data[:, icell, :], target_pressure)

        fields[name] = level

    fields["pressure_hPa"] = np.full(n_cells, pressure_hpa, dtype=np.float64)
    return centres, fields


def regular_lonlat_points():
    lons = np.linspace(lon_min, lon_max, lon_dim)
    lats = np.linspace(lat_min, lat_max, lat_dim)
    lon2d, lat2d = np.meshgrid(lons, lats)
    points = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    return points


def periodic_source_points(centres, values):
    """Duplicate source points at +/-360 deg to reduce dateline edge artefacts."""
    pts = np.vstack((centres, centres + np.array([360.0, 0.0]), centres - np.array([360.0, 0.0])))
    vals = np.concatenate((values, values, values), axis=0)
    return pts, vals


def horizontal_interpolate_field(centres, values, target_points, method="linear", fill_nearest=True):
    values = np.asarray(values, dtype=np.float64)

    if values.ndim == 1:
        ok = np.isfinite(values) & np.isfinite(centres[:, 0]) & np.isfinite(centres[:, 1])
        out = np.full(target_points.shape[0], np.nan, dtype=np.float64)
        if np.count_nonzero(ok) >= 3:
            pts, vals = periodic_source_points(centres[ok], values[ok])
            out = griddata(pts, vals, target_points, method=method)
            if fill_nearest and np.any(~np.isfinite(out)):
                near = griddata(pts, vals, target_points[~np.isfinite(out)], method="nearest")
                out[~np.isfinite(out)] = near
        elif np.count_nonzero(ok) > 0:
            pts, vals = periodic_source_points(centres[ok], values[ok])
            out = griddata(pts, vals, target_points, method="nearest")
        return out

    # Multi-component field.
    out = np.full((target_points.shape[0], values.shape[1]), np.nan, dtype=np.float64)
    for c in range(values.shape[1]):
        out[:, c] = horizontal_interpolate_field(centres, values[:, c], target_points, method, fill_nearest)
    return out


def make_pressure_level_image(vtp_series, pressure_hpa, horizontal_method="linear", fill_nearest=True):
    centres, fields = make_native_pressure_level_fields(vtp_series, pressure_hpa)
    target_points = regular_lonlat_points()

    image = vtk.vtkImageData()
    image.SetDimensions(lon_dim, lat_dim, 1)
    dx = (lon_max - lon_min) / max(lon_dim - 1, 1)
    dy = (lat_max - lat_min) / max(lat_dim - 1, 1)
    image.SetSpacing(dx, dy, 1.0)
    image.SetOrigin(lon_min, lat_min, float(pressure_hpa))

    for name, values in fields.items():
        grid_values = horizontal_interpolate_field(
            centres, values, target_points, method=horizontal_method, fill_nearest=fill_nearest
        )
        arr = numpy_to_vtk(grid_values, deep=True)
        arr.SetName(name)
        image.GetPointData().AddArray(arr)

    return image


def average_output_name(pressure_hpa):
    ptag = f"p{pressure_hpa:g}hPa".replace(".", "p")
    return f"{run}_{t1:04d}_{t2:04d}_mean_{ptag}.vti"


def accumulate_running_mean(mean_arrays, count_arrays, image):
    """Update per-array running means, ignoring NaNs point by point."""
    point_data = image.GetPointData()

    for i in range(point_data.GetNumberOfArrays()):
        name = point_data.GetArrayName(i)
        values = vtk_to_numpy(point_data.GetArray(i)).astype(np.float64)
        finite = np.isfinite(values)

        if name not in mean_arrays:
            mean_arrays[name] = np.zeros_like(values, dtype=np.float64)
            count_arrays[name] = np.zeros_like(values, dtype=np.int64)

        counts_old = count_arrays[name][finite]
        counts_new = counts_old + 1
        mean_arrays[name][finite] += (values[finite] - mean_arrays[name][finite]) / counts_new
        count_arrays[name][finite] = counts_new


def build_mean_image(mean_arrays, count_arrays, pressure_hpa):
    image = vtk.vtkImageData()
    image.SetDimensions(lon_dim, lat_dim, 1)
    dx = (lon_max - lon_min) / max(lon_dim - 1, 1)
    dy = (lat_max - lat_min) / max(lat_dim - 1, 1)
    image.SetSpacing(dx, dy, 1.0)
    image.SetOrigin(lon_min, lat_min, float(pressure_hpa))

    for name, mean in mean_arrays.items():
        data = mean.copy()
        data[count_arrays[name] == 0] = np.nan
        arr = numpy_to_vtk(data, deep=True)
        arr.SetName(name)
        image.GetPointData().AddArray(arr)
    return image


def write_image(image, filename):
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(image)
    ok = writer.Write()
    if ok != 1:
        raise RuntimeError(f"Failed to write {filename}")


class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


DOC = textwrap.dedent("""\
  Creates one time-averaged *.vti file over the requested time range.

  Each time index is first interpolated directly from native uneven pressure
  columns to a fixed pressure level, then horizontally interpolated to a regular
  lon-lat image and accumulated into a running average.

  Physical pressure is computed from the two data fields:

      P = (P/Ps) * Ps

  The default pressure level is 500 hPa.
""")


def parse_args():
    p = argparse.ArgumentParser(
        prog="lonlat_to_pressure_level_mean_native_griddata.py",
        description="Interpolate model output directly to a specified pressure level, average over time, and write one 2D lon-lat VTI file.",
        epilog=DOC,
        formatter_class=_Fmt,
    )

    p.add_argument("run")
    p.add_argument("Jmin", type=int)
    p.add_argument("Jmax", type=int)
    p.add_argument("nz", type=int)
    p.add_argument("t1", type=int)
    p.add_argument("t2", type=int, nargs="?", help="last time index for averaging; defaults to t1")
    p.add_argument("--pressure-hpa", type=float, default=500.0,
                   help="pressure level for interpolation, in hPa")
    p.add_argument("--horizontal-method", choices=["linear", "nearest"], default="linear",
                   help="horizontal interpolation method for native cell centres to regular lon-lat grid")
    p.add_argument("--no-fill-nearest", action="store_true",
                   help="leave linear-interpolation holes as NaN instead of filling with nearest neighbour")

    if len(sys.argv) == 1:
        print()
        p.print_help()
        sys.exit(0)

    a = p.parse_args()
    if a.t2 is None:
        a.t2 = a.t1

    if rank == 0:
        for k, v in sorted(vars(a).items()):
            print(f"{k:>18}: {v}", flush=True)

    return a


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

args = parse_args()

run = args.run
Jmin = args.Jmin
Jmax = args.Jmax
nz = args.nz
t1 = args.t1
t2 = args.t2

lat_dim = int(np.sqrt((10 * 4**Jmax + 2) / 2))
lon_dim = 2 * lat_dim

lon_min, lon_max = -180.0, 180.0
lat_min, lat_max = -90.0, 90.0

with suppress(OSError):
    os.remove(f"{run}/.DS_Store")

if rank == 0:
    print(f"\nDirectly interpolating native columns to {args.pressure_hpa:g} hPa", flush=True)
    print(f"Output grid: {lon_dim} x {lat_dim} lon-lat", flush=True)
    print(f"Horizontal interpolation: {args.horizontal_method}", flush=True)
    print(f"Averaging file indices {t1} to {t2}; only final mean VTI will be written\n", flush=True)

mean_arrays = {}
count_arrays = {}

for t in range(t1, t2 + 1):
    if rank == 0:
        print(f"Processing file index {t}", flush=True)

    failures = run_xyz_layers(run, Jmin, Jmax, nz, t, t, True)
    if failures and rank == 0:
        print("Failures:", failures, flush=True)

    comm.Barrier()

    if rank != 0:
        continue

    vtp_series = [f"{run}_tri_lonlat_{z:03d}_{t:04d}.vtp" for z in range(1, nz + 1)]
    pressure_level_image = make_pressure_level_image(
        vtp_series,
        args.pressure_hpa,
        horizontal_method=args.horizontal_method,
        fill_nearest=not args.no_fill_nearest,
    )
    accumulate_running_mean(mean_arrays, count_arrays, pressure_level_image)

    delete_files(f"{run}_tri*{t:04d}.vtp")

if rank == 0:
    if not mean_arrays:
        raise RuntimeError("No data were accumulated; check the requested time range and input files.")
    mean_image = build_mean_image(mean_arrays, count_arrays, args.pressure_hpa)
    filename = average_output_name(args.pressure_hpa)
    write_image(mean_image, filename)
    print(f"\nWrote time mean over {t1}..{t2}: {filename}", flush=True)

comm.Barrier()
