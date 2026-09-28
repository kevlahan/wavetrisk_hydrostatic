# Converts xyz data in a vtk file to lonlat data and saves in a new file
# (comutes a longitude-latitude projection of spherical shell data)
#
# Usage: python xyz2lonlat.py input_vtk_file output_vtk_file
#
# Weiguang Guan     (SHARCNET)            2024-04-16
# Nicholas Kevlahan (McMaster University) 2025-10-16
import argparse, os, glob, subprocess, sys, tarfile, textwrap
import hashlib
import numpy as np
import vtk
from contextlib import suppress
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk, numpy_to_vtkIdTypeArray
from pathlib import Path

verbose = False

class _Fmt(argparse.ArgumentDefaultsHelpFormatter,
           argparse.RawDescriptionHelpFormatter):
    """Show defaults and preserve newlines/indentation in epilog."""
    pass

def parse_args():
    epilog = textwrap.dedent("""\
      Converts vertical layers of <run>_tri_<t>.vtk.tgz to lon/lat .vtp.

      Examples
        # Single z,t
        python xyz2lonlat.py SimpleJ5Z30 5 5 1 1 365 365 y

        # Range of z with one time
        python xyz2lonlat.py SimpleJ5Z30 5 5 1 30 365 365 y

        # Range of z and t
        python xyz2lonlat.py SimpleJ5Z30 5 5 1 30 361 365 y

      Notes
        • Tarball expected in the working directory:
            <run>_tri_<t:04d>.vtk.tgz
        • The legacy Delaunay y/n flag is accepted; both use seam/pole splitting.
        • Combine with the MPI layer driver when parallelizing.
    """)
    p = argparse.ArgumentParser(
        prog="xyz2lonlat.py",
        description="Convert one or more layers/times from VTK to lon/lat VTP.",
        formatter_class=_Fmt,
        epilog=epilog,
    )
    # Positionals (kept to remain back-compatible with your current callers)
    p.add_argument("run", help="Run prefix (e.g., SimpleJ5Z30)")
    p.add_argument("Jmin", type=int, help="Minimum level (inclusive)")
    p.add_argument("Jmax", type=int, help="Maximum level (inclusive)")
    p.add_argument("z1",   type=int, help="First layer index")
    p.add_argument("z2",   type=int, help="Last  layer index (>= z1)")
    p.add_argument("t1",   type=int, help="First time index (used in file names)")
    p.add_argument("t2",   type=int, help="Last  time index (>= t1)")
    p.add_argument("Delaunay", choices=["y","n","Y","N"],
                   help="Legacy compatibility flag; both values use seam/pole splitting")

    if len(sys.argv) == 1:  # no args → show full help + examples
        print()
        p.print_help(); sys.exit(0)

    a = p.parse_args()
    
    if verbose:
        print(f"Converting layers for run={a.run}, J=[{a.Jmin},{a.Jmax}], "
              f"z2={a.z2}, z2={a.z2}, t1={a.t1}, t2={a.t2}, Delaunay={a.Delaunay}")

    # Unpack to variable names
    run          = a.run
    Jmin         = a.Jmin
    Jmax         = a.Jmax 
    z1           = a.z1
    z2           = a.z2
    t1           = a.t1
    t2           = a.t2
    Delaunay     = a.Delaunay

    return run, Jmin, Jmax, z1, z2, t1, t2, Delaunay


def _unwrap_longitudes(lon):
    """Put a local edge/triangle on its shortest continuous longitude interval."""
    wrapped = np.mod(lon, 360.0)
    ordered = np.sort(wrapped)
    gaps = np.diff(np.r_[ordered, ordered[0] + 360.0])
    start = ordered[(int(np.argmax(gaps)) + 1) % len(ordered)]
    return np.where(wrapped < start - 1e-12, wrapped + 360.0, wrapped)


def _clip_longitude(poly, bound, keep_greater):
    result = []
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        inside_a = a[0] >= bound if keep_greater else a[0] <= bound
        inside_b = b[0] >= bound if keep_greater else b[0] <= bound
        if inside_a:
            result.append(a.copy())
        if inside_a != inside_b:
            q = a + ((bound - a[0]) / (b[0] - a[0])) * (b - a)
            q[0] = bound
            result.append(q)
    return np.asarray(result, dtype=float).reshape(-1, 2)


_projection_cache = None  # One geometry per worker; fields are always refreshed.


def project_spherical_triangles(data):
    """Project a closed triangular shell without moving nearby physical vertices.

    Original cell values are copied to all child triangles. Longitude is
    periodic; seam-crossing faces are clipped into both sides of the map.
    Faces containing a pole are split at that pole, whose longitude is
    represented by a boundary segment rather than an arbitrary single point.
    Edges between nonpolar vertices remain straight in the lon/lat plane.
    This is a geometric projection, not a conservative spherical-area remap.
    """
    global _projection_cache
    xyz = vtk_to_numpy(data.GetPoints().GetData()).astype(np.float64)
    cell_array = data.GetPolys() if isinstance(data, vtk.vtkPolyData) else data.GetCells()
    offsets = vtk_to_numpy(cell_array.GetOffsetsArray())
    if not np.all(np.diff(offsets) == 3):
        raise ValueError("Expected triangular surface cells")
    if isinstance(data, vtk.vtkPolyData) and data.GetNumberOfCells() != cell_array.GetNumberOfCells():
        raise ValueError("Expected only triangular polygons")
    if isinstance(data, vtk.vtkUnstructuredGrid) and not np.all(vtk_to_numpy(data.GetCellTypesArray()) == vtk.VTK_TRIANGLE):
        raise ValueError("Expected only triangular surface cells")
    triangles = vtk_to_numpy(cell_array.GetConnectivityArray()).reshape(-1, 3)
    digest = hashlib.blake2b(digest_size=20)
    digest.update(np.ascontiguousarray(xyz).view(np.uint8))
    digest.update(np.ascontiguousarray(triangles, dtype=np.int64).view(np.uint8))
    key = (xyz.shape, triangles.shape, digest.digest())
    source_level = data.GetCellData().GetArray("Level")
    if source_level is not None:
        lv = vtk_to_numpy(source_level)
        if not np.isfinite(lv).all() or not np.allclose(lv, np.rint(lv), rtol=0, atol=1e-6):
            raise ValueError("Source cell Level must contain integers")

    def attach_fields(geometry, parents):
        out = vtk.vtkPolyData()
        out.CopyStructure(geometry)
        cell_data = data.GetCellData()
        if all(cell_data.GetArray(i) is not None for i in range(cell_data.GetNumberOfArrays())):
            for i in range(cell_data.GetNumberOfArrays()):
                source = cell_data.GetArray(i)
                values = vtk_to_numpy(source)[parents]
                if source.GetName() == "Level":
                    values = np.rint(values).astype(np.int32)
                array = numpy_to_vtk(values, deep=True)
                array.SetName(source.GetName() or f"Array_{i}")
                out.GetCellData().AddArray(array)
        else:
            out.GetCellData().CopyAllocate(cell_data, len(parents))
            for child, parent in enumerate(parents):
                out.GetCellData().CopyData(cell_data, int(parent), child)
        out.GetFieldData().ShallowCopy(data.GetFieldData())
        return out

    if _projection_cache is not None and _projection_cache[0] == key:
        return attach_fields(_projection_cache[1], _projection_cache[2])
    radii = np.linalg.norm(xyz, axis=1)
    if not np.isfinite(xyz).all() or np.any(radii <= 0):
        raise ValueError("Expected finite nonzero coordinates on a spherical shell")
    unit = xyz / radii[:, None]
    lon = np.degrees(np.arctan2(unit[:, 1], unit[:, 0]))
    lat = np.degrees(np.arctan2(unit[:, 2], np.linalg.norm(unit[:, :2], axis=1)))
    xy = np.column_stack((lon, lat))
    faces = unit[triangles]
    determinant = np.einsum('ij,ij->i', faces[:, 0], np.cross(faces[:, 1], faces[:, 2]))
    if np.any(np.abs(determinant) < 1e-14):
        raise ValueError("Degenerate spherical triangle")
    polar = np.linalg.norm(unit[:, :2], axis=1) < 1e-10
    ordinary = (np.ptp(lon[triangles], axis=1) < 180.0 - 1e-10) & ~polar[triangles].any(axis=1)
    # Most cells need only a coordinate transform. Only seam/pole cells enter Python loops.
    keys = np.round(xy, 10)
    _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    base_xy = xy[first]
    point_ids = {tuple(keys[i]): j for j, i in enumerate(first)}
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.column_stack((base_xy, np.zeros(len(base_xy)))), deep=True))
    ordinary_ids = np.flatnonzero(ordinary)
    regular = inverse[triangles[ordinary_ids]].copy()
    pxy = base_xy[regular]
    u, v = pxy[:, 1] - pxy[:, 0], pxy[:, 2] - pxy[:, 0]
    reverse = u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0] < 0
    regular[reverse] = regular[reverse][:, [0, 2, 1]]
    cells = vtk.vtkCellArray()  # Boundary triangles only.
    parents = []

    def insert_polygon(poly, parent):
        # Remove consecutive duplicates introduced by clipping at an endpoint.
        clean = []
        for point in poly:
            if not clean or np.linalg.norm(point - clean[-1]) > 1e-10:
                clean.append(point)
        if len(clean) > 1 and np.linalg.norm(clean[0] - clean[-1]) < 1e-10:
            clean.pop()
        if len(clean) < 3:
            return
        poly = np.asarray(clean)
        area = np.sum(poly[:, 0] * np.roll(poly[:, 1], -1) -
                      poly[:, 1] * np.roll(poly[:, 0], -1)) / 2
        if abs(area) < 1e-12:
            return
        if area < 0:
            poly = poly[::-1]
        ids = []
        for x, y in poly:
            # Only remove floating-point roundoff at constructed boundaries.
            if abs(abs(x) - 180.0) < 1e-10:
                x = np.copysign(180.0, x)
            if abs(abs(y) - 90.0) < 1e-10:
                y = np.copysign(90.0, y)
            # Weld shared edges within numerical tolerance; keep +/-180 distinct.
            key = (round(float(x), 10), round(float(y), 10))
            if key not in point_ids:
                point_ids[key] = points.InsertNextPoint(float(x), float(y), 0.0)
            ids.append(point_ids[key])
        for j in range(1, len(ids) - 1):
            tri = [ids[0], ids[j], ids[j+1]]
            p = np.array([points.GetPoint(i)[:2] for i in tri])
            u, v = p[1] - p[0], p[2] - p[0]
            if abs(u[0]*v[1] - u[1]*v[0]) < 1e-12:
                continue
            cells.InsertNextCell(3, tri)
            parents.append(parent)

    def map_face(face, parent):
        polar = np.flatnonzero(np.linalg.norm(face[:, :2], axis=1) < 1e-10)
        if polar.size:
            k = int(polar[0])
            others = [i for i in range(3) if i != k]
            edge = face[others]
            lon = _unwrap_longitudes(np.degrees(np.arctan2(edge[:, 1], edge[:, 0])))
            lat = np.degrees(np.arctan2(edge[:, 2], np.linalg.norm(edge[:, :2], axis=1)))
            pole_lat = float(np.copysign(90.0, face[k, 2]))
            poly = np.array([[lon[0], lat[0]], [lon[1], lat[1]],
                             [lon[1], pole_lat], [lon[0], pole_lat]])
        else:
            lon = _unwrap_longitudes(np.degrees(np.arctan2(face[:, 1], face[:, 0])))
            lat = np.degrees(np.arctan2(face[:, 2], np.linalg.norm(face[:, :2], axis=1)))
            poly = np.column_stack((lon, lat))
        # Every periodic copy which intersects the canonical map window.
        lower = int(np.ceil((-180.0 - poly[:, 0].max()) / 360.0))
        upper = int(np.floor((180.0 - poly[:, 0].min()) / 360.0))
        for shift in range(lower, upper + 1):
            clipped = poly.copy()
            clipped[:, 0] += 360.0 * shift
            clipped = _clip_longitude(clipped, -180.0, True)
            if len(clipped) >= 3:
                clipped = _clip_longitude(clipped, 180.0, False)
            if len(clipped) >= 3:
                insert_polygon(clipped, parent)

    for cid in np.flatnonzero(~ordinary):
        face = faces[cid]
        if np.any(np.linalg.norm(face[:, :2], axis=1) < 1e-10):
            map_face(face, int(cid))
            continue
        split = False
        for pole in (np.array([0., 0., 1.]), np.array([0., 0., -1.])):
            weights = np.linalg.solve(face.T, pole)
            if np.all(weights >= -1e-12):
                for i in range(3):
                    subface = np.array([face[i], face[(i+1) % 3], pole])
                    if abs(np.linalg.det(subface)) > 1e-14:
                        map_face(subface, int(cid))
                split = True
                break
        if not split:
            map_face(face, int(cid))
    boundary = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, 3)
    connectivity = np.concatenate((regular, boundary)).astype(np.int64)
    parents = np.concatenate((ordinary_ids, np.asarray(parents, dtype=np.int64)))
    # Remove unused original pole vertices, which would otherwise have no pressure support.
    used, new_ids = np.unique(connectivity, return_inverse=True)
    coords = vtk_to_numpy(points.GetData())[used]
    points.SetData(numpy_to_vtk(coords, deep=True))
    connectivity = new_ids.reshape(-1, 3).astype(np.int64)
    cells = vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, connectivity.size+1, 3, dtype=np.int64), deep=True),
                  numpy_to_vtkIdTypeArray(connectivity.ravel(), deep=True))
    out = vtk.vtkPolyData()
    out.SetPoints(points)
    out.SetPolys(cells)
    _projection_cache = (key, out, parents)
    return attach_fields(out, parents)


def convert_layer(run, z, t):
    """Convert an already extracted layer; callable by each MPI worker."""
    vtk_file = f"{run}_tri_{z:03d}_{t:04d}.vtk"
    vtp_file = f"{run}_tri_lonlat_{z:03d}_{t:04d}.vtp"
    reader = vtk.vtkDataSetReader()
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllFieldsOn()
    if not Path(vtk_file).is_file():
        raise FileNotFoundError(vtk_file)
    reader.SetFileName(vtk_file)
    reader.Update()
    data = reader.GetOutput()
    if data.GetNumberOfCells() == 0:
        raise RuntimeError(f"No cells read from {vtk_file}")
    output = project_spherical_triangles(data)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(vtp_file)
    writer.SetInputData(output)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write {vtp_file}")
    os.remove(vtk_file)
    return vtp_file



if __name__ == "__main__":
    (run, Jmin, Jmax, z1, z2, t1, t2, Delaunay) = parse_args()

    # Grid dimensions (same number of rectangular cells as lozenge cells on the sphere) 
    lat_dim = int (np.sqrt((10*4**Jmax + 2)/2))
    lon_dim = 2*lat_dim

    dtheta_min = 180/lat_dim
    dtheta_max = dtheta_min * 2**(Jmax - Jmin)

    wdir = str(Path.cwd().resolve())

    for t in range (t1, t2+1):
        for z in range (z1, z2+1):
            # Load the input vtk file
            vtk_file = f"{run}_tri_{z:03d}_{t:04d}.vtk"
            vtp_file = f"{run}_tri_lonlat_{z:03d}_{t:04d}.vtp"

            # Uncompress *.vtk.tgz if necessary
            p = Path(vtk_file).expanduser()
            if not p.exists():
                tar_file = f"{run}_tri_{t:04d}.vtk.tgz"

                p = Path(tar_file).expanduser()
                if not p.exists():
                    sys.exit(f"ERROR: file not found: {p}")   # exits with code 1
                    if verbose: print ("Uncompressing file ", file)

                with tarfile.open(tar_file, "r:gz") as tar:
                    tar.extractall(wdir, filter="data")

                # Remove surface data *.vtk file
                surface_vtk = f"{run}_tri_000_{t:04d}.vtk"
                p = Path(surface_vtk).expanduser()
                if p.exists():
                    os.remove(surface_vtk)

            convert_layer(run, z, t)
