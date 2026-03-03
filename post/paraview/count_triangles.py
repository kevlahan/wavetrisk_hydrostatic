#!/usr/bin/env python3
"""
Count triangle cells in legacy VTK (BINARY) POLYDATA written by Fortran
(unformatted stream, BIG_ENDIAN) inside *.vtk.tgz archives.

For each archive like:
  SimpleJ5J7Z30_tri_0365.vtk.tgz
it extracts only layer 1:
  SimpleJ5J7Z30_tri_001_0365.vtk
counts triangles from the POLYGONS connectivity, deletes extracted data,
and writes a CSV with columns: seq, n_triangles.
"""

import argparse
import csv
import re
import struct
import tarfile
import tempfile
from pathlib import Path


# Matches files ending with _####.vtk.tgz (e.g., ..._0365.vtk.tgz)
SAVE_RE = re.compile(r"_(\d{4})\.vtk\.tgz$")


def find_layer1_member(tar: tarfile.TarFile, save_num: int) -> tarfile.TarInfo:
    """
    Find the member corresponding to layer 1 for this save:
      *_001_####.vtk
    """
    suffix = f"_001_{save_num:04d}.vtk"
    candidates = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(suffix)]
    if not candidates:
        raise FileNotFoundError(f"No member ending with {suffix} found in archive {tar.name}")
    # If multiple, pick the shortest path name
    return min(candidates, key=lambda m: len(m.name))


def count_triangles_binary_polydata(vtk_path: Path) -> int:
    """
    Count triangles in a legacy VTK POLYDATA file produced by Fortran with:
      form='unformatted', access='stream', convert='BIG_ENDIAN'

    We do NOT assume the header is readable line-by-line as UTF-8/ASCII,
    because unformatted stream writes raw bytes.
    Instead:
      - read entire file bytes
      - locate the ASCII token 'POLYGONS'
      - parse the two ASCII integers on that same line: ncell, size
      - then read ncell polygon records from binary (big-endian int32):
            nvert, id1, id2, ..., id_nvert
      - count how many have nvert==3
    """
    data = vtk_path.read_bytes()

    key = b"POLYGONS"
    idx = data.find(key)
    if idx == -1:
        raise ValueError(f"{vtk_path.name}: POLYGONS section not found.")

    # Find end of the POLYGONS header line (ASCII newline)
    line_end = data.find(b"\n", idx)
    if line_end == -1:
        raise ValueError(f"{vtk_path.name}: malformed POLYGONS line (no newline).")

    # Decode just that header line safely as ASCII
    line = data[idx:line_end].decode("ascii", errors="strict")
    parts = line.split()
    if len(parts) < 3 or parts[0].upper() != "POLYGONS":
        raise ValueError(f"{vtk_path.name}: unexpected POLYGONS header: {line!r}")

    n_polys = int(parts[1])
    poly_size = int(parts[2])

    # Binary connectivity begins immediately after the newline
    offset = line_end + 1

    be_i4 = struct.Struct(">i")  # big-endian 32-bit int

    n_tri = 0
    ints_read = 0

    for _ in range(n_polys):
        if offset + 4 > len(data):
            raise EOFError(f"{vtk_path.name}: EOF while reading polygon nvert.")
        nvert = be_i4.unpack_from(data, offset)[0]
        offset += 4
        ints_read += 1

        if nvert == 3:
            n_tri += 1

        skip = 4 * nvert
        if offset + skip > len(data):
            raise EOFError(f"{vtk_path.name}: EOF while skipping polygon indices.")
        offset += skip
        ints_read += nvert

    # Sanity check (VTK meaning of size is total number of ints in connectivity)
    if ints_read != poly_size:
        print(f"Warning: {vtk_path.name}: POLYGONS size={poly_size} but parsed ints={ints_read}")

    return n_tri


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Count triangles from layer-1 legacy VTK BINARY POLYDATA inside *.vtk.tgz archives and write CSV."
    )
    ap.add_argument("input_dir", type=Path, help="Directory containing *.vtk.tgz files")
    ap.add_argument("-o", "--output", type=Path, default=Path("triangle_counts.csv"),
                    help="Output CSV filename (default: triangle_counts.csv)")
    ap.add_argument("--glob", default="*.vtk.tgz", help="Glob pattern for archives (default: *.vtk.tgz)")
    args = ap.parse_args()

    archives = sorted(args.input_dir.glob(args.glob))

    # Map save number -> archive path
    save_to_archive: dict[int, Path] = {}
    for p in archives:
        m = SAVE_RE.search(p.name)
        if m:
            save_to_archive[int(m.group(1))] = p

    if not save_to_archive:
        raise SystemExit(f"No archives matching *_####.vtk.tgz found in {args.input_dir}")

    results: list[tuple[int, int]] = []

    for save_num in sorted(save_to_archive.keys()):
        tgz_path = save_to_archive[save_num]
        print(f"Processing save {save_num:04d}: {tgz_path.name}")

        with tarfile.open(tgz_path, "r:gz") as tar:
            member = find_layer1_member(tar, save_num)

            # Extract only this member to a temp directory; temp directory is deleted automatically
            with tempfile.TemporaryDirectory(prefix=f"vtk_extract_{save_num:04d}_") as tmpd:
                tmpd_path = Path(tmpd)

                # Avoid Python 3.14+ tarfile extraction behavior changes by using a filter
                try:
                    tar.extract(member, path=tmpd_path, filter="data")
                except TypeError:
                    # Older Python: no 'filter' argument
                    tar.extract(member, path=tmpd_path)

                vtk_path = (tmpd_path / member.name).resolve()
                n_tri = count_triangles_binary_polydata(vtk_path)
                results.append((save_num, n_tri))

    # Write CSV
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["seq", "n_triangles"])
        for seq, ntri in results:
            w.writerow([seq, ntri])

    print(f"Wrote {len(results)} rows to: {args.output}")


if __name__ == "__main__":
    main()
