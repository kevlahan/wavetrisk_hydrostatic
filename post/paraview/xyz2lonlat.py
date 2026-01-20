# Converts xyz data in a vtk file to lonlat data and saves in a new file
# (comutes a longitude-latitude projection of spherical shell data)
#
# Usage: python xyz2lonlat.py input_vtk_file output_vtk_file
#
# Weiguang Guan     (SHARCNET)            2024-04-16
# Nicholas Kevlahan (McMaster University) 2025-10-16
import argparse, os, glob, subprocess, sys, tarfile, textwrap
import numpy as np
import vtk
from contextlib import suppress
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
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
        • use_flag: 'y' to apply Delaunay2D (gap fill), 'n' to skip.
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
                   help="y => Use Delaunay2D gap-fill")

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

    return run, Jmin, Jmax, z1, z2, t1, t2, t1, t2, Delaunay


if __name__ == "__main__":
    (run, Jmin, Jmax, z1, z2, t1, t2, t1, t2, Delaunay) = parse_args()

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

            vtkreader = vtk.vtkDataSetReader()
            vtkreader.ReadAllScalarsOn()
            vtkreader.SetFileName(vtk_file)
            vtkreader.Update()

            data = vtkreader.GetOutput()

            # Get point data and convert them from (x,y,z) to (lon, lat, 0)
            points     = data.GetPoints()
            coords     = vtk_to_numpy(points.GetData())
            num_points = coords.shape[0]

            # Compute radius of sphere
            R = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)

            # Conversion from x,y,z to lon, lat, 0
            coords[:,0] = np.degrees(np.arctan2(coords[:,1], coords[:,0])) # longitude
            coords[:,1] = np.degrees(np.arcsin(coords[:,2] / R))           # latitude
            coords[:,2] = 0

            points.SetData(numpy_to_vtk(coords))

            # Get cell array
            if isinstance(data, vtk.vtkPolyData):
                cells = data.GetPolys()
            elif isinstance(data, vtk.vtkUnstructuredGrid):
                cells = data.GetCells()

            cellformation = vtk_to_numpy(cells.GetData())
            num_cells     = data.GetNumberOfCells()

            # Create point to cell mapping and an array of start indices of cells
            point2cells = [ [] for _ in range(num_points) ] # a list of lists
            startIDs    = [0] * num_cells # start index in cell data "cellformation"

            startID = 0
            for cell in range(num_cells):              # loop through cells
                startIDs[cell] = startID
                size           = cellformation[startID] # number of vertices
                for i in range(size):
                    pid = cellformation[startID+1+i]
                    point2cells[pid].append(cell)

                startID = startID + 1 + size

            # Loop over cells, check if it intersects with -180/180
            # separation line. If yes, move its vertices on one side horizontally
            # so that their lon coordinates are -180. Which side's vertices to be moved 
            # depends on numbers of vertices of that cell on both sides.
            points_on_sep = set()
            for cell in range(num_cells):         # loop through cells
                startID  = startIDs[cell]
                size     = cellformation[startID] # number of vertices                
                pids_pos = []
                pids_neg = []
                for i in range(size):
                    pid = cellformation[startID+1+i]
                    if (  coords[pid,0]<-90): # on negative side
                        pids_neg.append(pid)
                    elif (coords[pid,0]> 90): # on positive side
                        pids_pos.append(pid)

                if (len(pids_pos)>len(pids_neg)):    # positive triangle
                    points_on_sep.update(pids_neg)
                else:                                # negative triangle
                    points_on_sep.update(pids_pos)
                    startIDs[cell] = -startIDs[cell] # sign used to indicate which side cell should be placed on

            # Split points on separation line
            last_pid = num_points
            for pid in points_on_sep:
                # Add a new point
                if (coords[pid,0]>0):
                    new_lon = -360 + coords[pid,0]
                else:
                    new_lon =  360 + coords[pid,0]

                points.InsertNextPoint(new_lon, coords[pid,1], coords[pid,2])
                coords = np.vstack((coords, [[new_lon, coords[pid,1], coords[pid,2]]]))

                # Let newly added point be vertices of some triangles
                for cid in point2cells[pid]: # loop over all cells sharing point pid
                    # Replace pid with last_pid for those cells sharing pid
                    startID = abs(startIDs[cid])
                    size = cellformation[startID]
                    if (startIDs[cid]>0 and coords[pid,0]<0) or (startIDs[cid]<0 and coords[pid,0]>0): 
                        loc = np.where(cellformation[startID+1:startID+1+size]==pid)
                        loc  = int(loc[0][0])
                        cells.ReplaceCellPointAtId(cid, loc, last_pid)
                        cellformation[startID+1+loc] = last_pid
                last_pid += 1

            points.SetData(numpy_to_vtk(coords))

            # Ensure no vertices are outside [-180,180] and fill gaps at poles
            num_points = coords.shape[0]
            for pid in range(num_points):
                if (np.abs(coords[pid,0])>180):
                    coords[pid,0] = np.sign(coords[pid,0]) * 180

            # Convert to point data
            cellToPoint = vtk.vtkCellDataToPointData()
            cellToPoint.SetInputData(data)
            cellToPoint.Update()
            output = cellToPoint.GetOutput()

            if Delaunay=='y':
                # Apply vtkDelaunay2D filter
                delaunay = vtk.vtkDelaunay2D()
                delaunay.SetInputData(output)
                delaunay.Update()
                output = delaunay.GetOutput()

                point2cell = vtk.vtkPointDataToCellData()
                point2cell.SetInputData(output)
                point2cell.Update()
                polydata = point2cell.GetOutput()

                # Shift points close to longitude boundaries to edges
                points = polydata.GetPoints()
                for i in range(polydata.GetNumberOfCells()):
                    cell   = polydata.GetCell(i)
                    pt_ids = cell.GetPointIds()
                    p      = np.array([points.GetPoint(pt_ids.GetId(0)), points.GetPoint(pt_ids.GetId(1)), points.GetPoint(pt_ids.GetId(2))])
                    dtheta = np.min([np.linalg.norm(p[0]-p[1]), np.linalg.norm(p[0]-p[2]), np.linalg.norm(p[1]-p[2])])
                    for j in range(3):
                        pid = pt_ids.GetId(j)
                        coord = list(points.GetPoint(pid))
                        if (np.abs(np.abs(coord[0])-180)<2*dtheta):  # modify to avoid gaps if necessary
                            coord[0] = np.sign(coord[0]) * 180
                            points.SetPoint(pid, coord)
                        if (np.abs(np.abs(coord[1])-90)<2*dtheta): # modify to avoid gaps if necessary
                            coord[1] = np.sign(coord[1]) * 90
                            points.SetPoint(pid, coord)
                points.Modified()
            else:
                point2cell = vtk.vtkPointDataToCellData()
                point2cell.SetInputData(output)
                point2cell.Update()
                polydata = point2cell.GetOutput()

            # Write out structured data
            writer = vtk.vtkXMLPolyDataWriter() 
            writer.SetFileName(vtp_file)
            writer.SetInputData(polydata)
            writer.Write()

            # Remove vtk file
            os.remove(vtk_file)

  
