#!/usr/bin/env python3
import argparse
import csv
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_vti(path: str):
    r = vtk.vtkXMLImageDataReader()
    r.SetFileName(path)
    r.Update()
    return r.GetOutput()


def point_dims(img):
    nx, ny, nz = img.GetDimensions()
    if nx <= 0 or ny <= 0 or nz <= 1:
        raise ValueError(f"Need at least 2 points in z for point-data z-integration. dims={img.GetDimensions()}")
    return nx, ny, nz


def latitude_values_point(img):
    """Point-coordinate latitude values y."""
    ox, oy, oz = img.GetOrigin()
    sx, sy, sz = img.GetSpacing()
    nx, ny, nz = point_dims(img)
    return oy + np.arange(ny) * sy


def z_integral_profile_y_point(img, vtk_arr):
    """
    PointData array shaped (nx, ny, nz, ncomp).
    Returns integrated profile vs y: shape (ny, ncomp).

    Uses trapezoidal rule in z with uniform spacing dz.
    If nx == 1, drops x. If nx > 1, averages over x after z-integration.
    """
    sx, sy, sz = img.GetSpacing()
    dz = float(sz)

    nx, ny, nz = point_dims(img)
    ncomp = vtk_arr.GetNumberOfComponents()

    a = vtk_to_numpy(vtk_arr)
    if ncomp == 1:
        a = a.reshape(-1, 1)

    # Reshape to structured. VTK ImageData uses Fortran-like ordering for reshape.
    a4 = a.reshape((nx, ny, nz, ncomp), order="F")

    # Trapezoidal weights along z
    w = np.ones(nz, dtype=np.float64)
    w[0] = 0.5
    w[-1] = 0.5

    # Integral over z -> (nx, ny, ncomp)
    Iz = (a4 * w[None, None, :, None]).sum(axis=2) * dz

    # Reduce x if needed
    if nx == 1:
        prof = Iz[0, :, :]          # (ny, ncomp)
    else:
        prof = Iz.mean(axis=0)      # (ny, ncomp)

    return prof


def main():
    ap = argparse.ArgumentParser(
        description="PointData-only: integrate all fields over z (trapz) and save CSV vs latitude (y)."
    )
    ap.add_argument("input_vti")
    ap.add_argument("output_csv")
    ap.add_argument("--lat-name", default="latitude", help="Latitude column name in CSV.")
    args = ap.parse_args()

    img = read_vti(args.input_vti)

    pd = img.GetPointData()
    if pd is None or pd.GetNumberOfArrays() == 0:
        raise RuntimeError("No PointData arrays found in this .vti.")

    nx, ny, nz = point_dims(img)
    lat = latitude_values_point(img)

    columns = {}
    for i in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(i)
        name = arr.GetName() or f"array_{i}"

        prof = z_integral_profile_y_point(img, arr)  # (ny, ncomp)
        ncomp = arr.GetNumberOfComponents()

        if prof.shape[0] != ny:
            raise RuntimeError(f"Latitude length mismatch for {name}: got {prof.shape[0]}, expected {ny}")

        if ncomp == 1:
            columns[name] = prof[:, 0]
        else:
            for j in range(ncomp):
                columns[f"{name}[{j}]"] = prof[:, j]

    # Write CSV
    with open(args.output_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        header = [args.lat_name] + list(columns.keys())
        w.writerow(header)
        for j in range(ny):
            w.writerow([float(lat[j])] + [float(columns[k][j]) for k in columns.keys()])

    print(f"Wrote: {args.output_csv}")


if __name__ == "__main__":
    main()
