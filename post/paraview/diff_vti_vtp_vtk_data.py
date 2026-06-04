#!/usr/bin/env pvpython

import argparse
import sys
import re
import os
import numpy as np

from vtkmodules.vtkIOXML import (
    vtkXMLImageDataReader,
    vtkXMLImageDataWriter,
    vtkXMLPolyDataReader,
    vtkXMLPolyDataWriter,
)
from vtkmodules.vtkIOLegacy import vtkDataSetReader, vtkDataSetWriter
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk


def dataset_type_from_extension(path):
    ext = os.path.splitext(path)[1].lower()

    if ext == ".vti":
        return "vti", vtkXMLImageDataReader, vtkXMLImageDataWriter
    if ext == ".vtp":
        return "vtp", vtkXMLPolyDataReader, vtkXMLPolyDataWriter
    if ext == ".vtk":
        return "vtk", vtkDataSetReader, vtkDataSetWriter

    raise ValueError(f"Unsupported extension {ext!r}; expected .vti, .vtp, or .vtk")


def read_dataset(path, Reader, dtype):
    r = Reader()
    r.SetFileName(path)

    # Important for legacy .vtk files: otherwise only active/default arrays
    # may be read.
    if dtype == "vtk":
        r.ReadAllScalarsOn()
        r.ReadAllVectorsOn()
        r.ReadAllTensorsOn()
        r.ReadAllNormalsOn()
        r.ReadAllColorScalarsOn()
        r.ReadAllFieldsOn()

    r.Update()
    return r.GetOutput()


def dataset_signature(data, dtype):
    if dtype == "vti":
        return (
            data.GetClassName(),
            tuple(data.GetExtent()),
            tuple(data.GetOrigin()),
            tuple(data.GetSpacing()),
            data.GetNumberOfPoints(),
            data.GetNumberOfCells(),
        )

    return (
        data.GetClassName(),
        data.GetNumberOfPoints(),
        data.GetNumberOfCells(),
    )


def list_arrays(dobj):
    return [
        (dobj.GetArray(i).GetName(), dobj.GetArray(i))
        for i in range(dobj.GetNumberOfArrays())
    ]


def canon(name):
    if name is None:
        return ""
    s = name.strip().lower()
    s = re.sub(r"[_\-\s/]+", "", s)
    s = re.sub(r"[^0-9a-z]+", "", s)
    return s


def build_map(dobj, label, strict_dupes=False):
    m = {}

    for name, arr in list_arrays(dobj):
        k = canon(name)

        if k in m:
            msg = (
                f"{label}: duplicate canonical key {k!r} for "
                f"{m[k][0]!r} and {name!r}. Keeping first."
            )
            if strict_dupes:
                raise RuntimeError(msg)
            print("WARNING:", msg, file=sys.stderr)
            continue

        m[k] = (name, arr)

    return m


def subtract_arrays(arrA, arrB, out_name):
    nA, nB = arrA.GetNumberOfTuples(), arrB.GetNumberOfTuples()
    cA, cB = arrA.GetNumberOfComponents(), arrB.GetNumberOfComponents()

    if nA != nB or cA != cB:
        raise ValueError(
            f"shape mismatch tuples/components: "
            f"A=({nA},{cA}) B=({nB},{cB})"
        )

    a = vtk_to_numpy(arrA)
    b = vtk_to_numpy(arrB)

    dt = np.result_type(a.dtype, b.dtype, np.float32)
    d = a.astype(dt, copy=False) - b.astype(dt, copy=False)

    out = numpy_to_vtk(d, deep=True)
    out.SetNumberOfComponents(cA)
    out.SetName(out_name)

    return out


def normalized_errors(a, b):
    na = vtk_to_numpy(a)
    nb = vtk_to_numpy(b)

    if na.shape != nb.shape:
        raise ValueError(f"shape mismatch {na.shape} vs {nb.shape}")

    da = na.ravel()
    db = nb.ravel()
    diff = da - db

    diff_f = diff.astype(np.float64, copy=False)
    da_f = da.astype(np.float64, copy=False)

    num_l2 = float(np.linalg.norm(diff_f))
    den_l2 = float(np.linalg.norm(da_f))

    num_inf = float(np.max(np.abs(diff_f))) if diff_f.size else 0.0
    den_inf = float(np.max(np.abs(da_f))) if da_f.size else 0.0

    def safe_div(num, den):
        if den == 0.0:
            return 0.0 if num == 0.0 else float("inf")
        return num / den

    return safe_div(num_l2, den_l2), safe_div(num_inf, den_inf)


def process(label, DA_A, DA_B, DA_out, args, err_fh):
    mapA = build_map(DA_A, label + " A", strict_dupes=args.strict_dupes)
    mapB = build_map(DA_B, label + " B", strict_dupes=args.strict_dupes)

    keys = sorted(set(mapA.keys()) & set(mapB.keys()))

    if args.list:
        print(f"{label} matched arrays ({len(keys)}):")
        for k in keys:
            nA, a = mapA[k]
            nB, b = mapB[k]
            print(
                f"  {nA!r}  <->  {nB!r}   "
                f"(key={k!r}, tuples={a.GetNumberOfTuples()}, "
                f"components={a.GetNumberOfComponents()})"
            )
        return 0

    ok = 0
    skip = 0

    for k in keys:
        nA, a = mapA[k]
        nB, b = mapB[k]

        try:
            DA_out.AddArray(subtract_arrays(a, b, nA))

            l2n, linfn = normalized_errors(a, b)
            if err_fh is not None:
                err_fh.write(
                    f"{label + ':' + nA:<40} "
                    f"{l2n:>14.8e} "
                    f"{linfn:>14.8e}\n"
                )

            ok += 1

        except Exception as e:
            skip += 1
            msg = f"{label}: skipping {nA!r} <-> {nB!r}: {e}"
            if args.strict:
                raise RuntimeError(msg) from e
            print("WARNING:", msg, file=sys.stderr)

    onlyA = sorted(set(mapA.keys()) - set(mapB.keys()))
    onlyB = sorted(set(mapB.keys()) - set(mapA.keys()))

    print(f"{label}: wrote {ok} arrays (skipped {skip}).")

    if onlyA:
        print(
            f"{label}: unmatched in A ({len(onlyA)}): "
            + ", ".join(mapA[k][0] for k in onlyA),
            file=sys.stderr,
        )

    if onlyB:
        print(
            f"{label}: unmatched in B ({len(onlyB)}): "
            + ", ".join(mapB[k][0] for k in onlyB),
            file=sys.stderr,
        )

    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Subtract matching arrays in A-B for .vtk, .vti, or .vtp files."
    )
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["point", "cell", "field", "all"], default="all")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--strict-dupes", action="store_true")
    ap.add_argument(
        "--ascii",
        action="store_true",
        help="For legacy .vtk output, write ASCII instead of binary.",
    )

    args = ap.parse_args()

    dtypeA, ReaderA, _ = dataset_type_from_extension(args.A)
    dtypeB, ReaderB, _ = dataset_type_from_extension(args.B)
    dtypeO, _, WriterO = dataset_type_from_extension(args.out)

    if dtypeA != dtypeB:
        raise RuntimeError(
            f"Input types differ: {args.A} is .{dtypeA}, {args.B} is .{dtypeB}"
        )

    if dtypeO != dtypeA:
        raise RuntimeError(
            f"Output extension must match input type: expected .{dtypeA}"
        )

    A = read_dataset(args.A, ReaderA, dtypeA)
    B = read_dataset(args.B, ReaderB, dtypeB)

    if args.mode in ("point", "cell", "all"):
        sigA = dataset_signature(A, dtypeA)
        sigB = dataset_signature(B, dtypeB)

        if sigA != sigB:
            msg = (
                "Dataset mismatch. Geometry/topology signatures differ.\n"
                f"A signature: {sigA}\n"
                f"B signature: {sigB}"
            )

            if args.strict:
                raise RuntimeError(msg)

            print("ERROR:", msg, file=sys.stderr)
            return 2

    output = A.NewInstance()
    output.ShallowCopy(A)

    # Write only difference arrays, not original arrays.
    output.GetPointData().Initialize()
    output.GetCellData().Initialize()
    output.GetFieldData().Initialize()

    base, _ = os.path.splitext(args.out)
    err_path = base + "_errors.txt"

    err_fh = None
    if not args.list:
        err_fh = open(err_path, "w", encoding="utf-8")
        err_fh.write("Field                                   rel l2 err    rel linf err\n")

    try:
        if args.mode in ("point", "all"):
            process(
                "PointData",
                A.GetPointData(),
                B.GetPointData(),
                output.GetPointData(),
                args,
                err_fh,
            )

        if args.mode in ("cell", "all"):
            process(
                "CellData",
                A.GetCellData(),
                B.GetCellData(),
                output.GetCellData(),
                args,
                err_fh,
            )

        if args.mode in ("field", "all"):
            process(
                "FieldData",
                A.GetFieldData(),
                B.GetFieldData(),
                output.GetFieldData(),
                args,
                err_fh,
            )

        if args.list:
            return 0

        w = WriterO()
        w.SetFileName(args.out)
        w.SetInputData(output)

        if dtypeO == "vtk":
            if args.ascii:
                w.SetFileTypeToASCII()
            else:
                w.SetFileTypeToBinary()

        if w.Write() != 1:
            raise RuntimeError(f"Failed to write {args.out}")

        print(f"Wrote: {args.out}")

        if err_fh is not None:
            err_fh.flush()
            err_fh.close()
            print(f"Wrote: {err_path}")

        return 0

    finally:
        if err_fh is not None and not err_fh.closed:
            err_fh.close()


if __name__ == "__main__":
    raise SystemExit(main())
