#!/usr/bin/env pvpython
import argparse
import sys
import re
import os
import numpy as np

from vtkmodules.vtkIOXML import vtkXMLImageDataReader, vtkXMLImageDataWriter
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk


def read_vti(path):
    r = vtkXMLImageDataReader()
    r.SetFileName(path)
    r.Update()
    return r.GetOutput()


def grid_signature(img):
    return (tuple(img.GetExtent()), tuple(img.GetOrigin()), tuple(img.GetSpacing()))


def list_arrays(dobj):
    return [(dobj.GetArray(i).GetName(), dobj.GetArray(i)) for i in range(dobj.GetNumberOfArrays())]


def canon(name: str) -> str:
    """
    Canonicalize an array name for matching:
    - lower case
    - remove underscores, spaces, hyphens, slashes
    - remove any non-alphanumeric characters
    """
    if name is None:
        return ""
    s = name.strip().lower()
    s = re.sub(r"[_\-\s/]+", "", s)        # drop common separators
    s = re.sub(r"[^0-9a-z]+", "", s)       # drop anything else non-alphanumeric
    return s


def build_map(dobj, label, strict_dupes=False):
    """
    Return dict: canonical_key -> (original_name, vtkArray)
    Warn on duplicates (after canonicalization).
    """
    m = {}
    for name, arr in list_arrays(dobj):
        k = canon(name)
        if k in m:
            msg = (f"{label}: duplicate canonical key {k!r} for "
                   f"{m[k][0]!r} and {name!r}. Keeping first.")
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
        raise ValueError(f"shape mismatch tuples/components: A=({nA},{cA}) B=({nB},{cB})")

    a = vtk_to_numpy(arrA)
    b = vtk_to_numpy(arrB)
    dt = np.result_type(a.dtype, b.dtype, np.float32)
    d = a.astype(dt, copy=False) - b.astype(dt, copy=False)

    out = numpy_to_vtk(d, deep=True)
    out.SetNumberOfComponents(cA)
    out.SetName(out_name)   # keep A's original field name
    return out

def normalized_errors(a, b):
    """
    Return (l2_norm, linf_norm) where:
      l2_norm  = ||a-b||_2 / ||a||_2
      linf_norm = ||a-b||_inf / ||a||_inf

    Arrays may be multi-component; treated as flattened vectors.
    If the denominator is 0, returns inf if numerator>0 else 0.
    """
    na = vtk_to_numpy(a)
    nb = vtk_to_numpy(b)

    if na.shape != nb.shape:
        raise ValueError(f"shape mismatch {na.shape} vs {nb.shape}")

    # Flatten including components
    da = na.ravel()
    db = nb.ravel()
    diff = da - db

    # Promote to float for norms even if inputs are ints
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["point", "cell", "field", "all"], default="all")
    ap.add_argument("--strict", action="store_true", help="Fail on any mismatch instead of skipping.")
    ap.add_argument("--list", action="store_true", help="List matched pairs and exit.")
    ap.add_argument("--strict-dupes", action="store_true",
                    help="Fail if canonicalization causes duplicate-name collisions within a dataset.")
    args = ap.parse_args()

    A = read_vti(args.A)
    B = read_vti(args.B)

    # Grid must match for point/cell differencing.
    if args.mode in ("point", "cell", "all"):
        if grid_signature(A) != grid_signature(B):
            msg = "Grid mismatch (extent/origin/spacing differ). Resample first."
            if args.strict:
                raise RuntimeError(msg)
            print("ERROR:", msg, file=sys.stderr)
            return 2

    # Output geometry only; wipe arrays so we write ONLY diffs
    out_img = vtkImageData()
    out_img.ShallowCopy(A)
    out_img.GetPointData().Initialize()
    out_img.GetCellData().Initialize()
    out_img.GetFieldData().Initialize()

    # ---- NEW: open an errors file unless --list ----
    base, _ = os.path.splitext(args.out)
    err_path = base + "_errors.txt"
    err_fh = None
    if not args.list:
        err_fh = open(err_path, "w", encoding="utf-8")
        err_fh.write("Field                             rel l2 err    rel linf err\n")

    def process(label, DA_A, DA_B, DA_out, err_fh):
        mapA = build_map(DA_A, label + " A", strict_dupes=args.strict_dupes)
        mapB = build_map(DA_B, label + " B", strict_dupes=args.strict_dupes)

        keys = sorted(set(mapA.keys()) & set(mapB.keys()))
        if args.list:
            print(f"{label} matched arrays ({len(keys)}):")
            for k in keys:
                nA, a = mapA[k]
                nB, b = mapB[k]
                print(f"  {nA!r}  <->  {nB!r}   (key={k!r}, tuples={a.GetNumberOfTuples()}, comps={a.GetNumberOfComponents()})")
            return 0

        ok = skip = 0
        for k in keys:
            nA, a = mapA[k]
            nB, b = mapB[k]
            try:
                # write diff array (existing behavior)
                DA_out.AddArray(subtract_arrays(a, b, nA))

                # ---- compute and record normalized errors ----
                l2n, linfn = normalized_errors(a, b)
                if err_fh is not None:
                    err_fh.write(
                        f"{nA:<32} "
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

        # Also report what did NOT match (helps catch typos)
        onlyA = sorted(set(mapA.keys()) - set(mapB.keys()))
        onlyB = sorted(set(mapB.keys()) - set(mapA.keys()))
        print(f"{label}: wrote {ok} arrays (skipped {skip}).")
        if onlyA:
            print(f"{label}: unmatched in A ({len(onlyA)}): " +
                  ", ".join(mapA[k][0] for k in onlyA), file=sys.stderr)
        if onlyB:
            print(f"{label}: unmatched in B ({len(onlyB)}): " +
                  ", ".join(mapB[k][0] for k in onlyB), file=sys.stderr)

        return 0

    try:
        if args.mode in ("point", "all"):
            process("PointData", A.GetPointData(), B.GetPointData(), out_img.GetPointData(), err_fh)

        if args.mode in ("cell", "all"):
            process("CellData", A.GetCellData(), B.GetCellData(), out_img.GetCellData(), err_fh)

        if args.mode in ("field", "all"):
            process("FieldData", A.GetFieldData(), B.GetFieldData(), out_img.GetFieldData(), err_fh)

        if args.list:
            return 0

        w = vtkXMLImageDataWriter()
        w.SetFileName(args.out)
        w.SetInputData(out_img)
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

