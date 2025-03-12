# Converts xyz data in a vtk file to lonlat data and saves in a new file
# (comutes a longitude-latitude projection of spherical shell data)
#
# Usage: python xyz2lonlat.py input_vtk_file output_vtk_file
#
# Weiguang Guan     (SHARCNET)            2024-04-16
# Nicholas Kevlahan (McMaster University) 2025-01-29
import os
import glob
import sys
import numpy as np
from utilities import *
import vtk
import subprocess
from contextlib import suppress
import tarfile
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

def safe_extract(tar, path=".", members=None):
    # Safely extract files, ensuring no files are extracted outside the target directory.

    for member in tar.getmembers():
        if not os.path.abspath(os.path.join(path, member.name)).startswith(os.path.abspath(path)):
            raise Exception(f"Unsafe extraction attempt: {member.name}")
    tar.extractall(path, members, filter="data")
    
def untar_files(t) :
    # Untars time t data
    
    file = run+'_tri_'+str(t).zfill(4)+".vtk.tgz"
    
    directory = os.getcwd()
    output_directory = directory
    
    tar_path = os.path.join(directory, file)
    try:
        with tarfile.open(tar_path, 'r:*') as tar:
            safe_extract(tar, path=output_directory)
    except tarfile.TarError as e:
        print(f"    Error extracting {file}: {e}")
    except Exception as e:
        print(f"    Security issue extracting {file}: {e}")
                    
# Main program
if (len(sys.argv)<7) :
    print("\nUsage: python xyz2lonlat.py run Jmin Jmax z1 z2 t1 t2 Delaunay\n")
    print("run      = file base name (without tri)")
    print("Jmin     = minimum level")
    print("Jmax     = maximum level")
    print("z1       = first z layer")
    print("z2       = last  z layer")
    print("t1       = first time")
    print("t2       = last time")
    print("Delaunay = 'y'/'n' (interpolate to Delaunay grid to remove gaps)\n")
    print("output is file_lonlat run_tri_zzz_tttt.vtk or or run_tri_zzz_tttt.vtp (Delaunay == 'y')\n")
    print("Example:")
    print("python3 xyz2lonlat.py SimpleJ5J1Z30 5 7 1 30 0 100 y")
    exit(0)

run       = sys.argv[1]
Jmin      = int(sys.argv[2])
Jmax      = int(sys.argv[3])
z1        = int(sys.argv[4])
z2        = int(sys.argv[5])
t1        = int(sys.argv[6])
t2        = int(sys.argv[7])
Delaunay  = sys.argv[8]

N        = int(np.sqrt(20*4**Jmax))
lat_dim  = int(N/2)
lon_dim  = 2*lat_dim

dtheta_min = 360.0/lon_dim
dtheta_max = dtheta_min * 2**(Jmax-Jmin)

straight_edges = False # fix straight edges on longitude boundaries

for t in range (t1, t2+1):
    untar_files(t)
    for z in range (z1, z2+1):
        # Load the input vtk file
        infile  = run+"_tri_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        outfile = run+"_tri_lonlat_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        #print("Transforming file "+infile+".vtk")
        vtkreader = vtk.vtkUnstructuredGridReader()
        vtkreader.ReadAllScalarsOn()
        vtkreader.SetFileName(infile+".vtk")
        vtkreader.Update()

        # Get the unstructed grid data
        ugrid = vtkreader.GetOutput()

        # Get point data and convert them from (x,y,z) to (lon, lat, 0)
        points     = ugrid.GetPoints()
        coords     = vtk_to_numpy(points.GetData())
        num_points = coords.shape[0]

        # Compute radius of sphere
        R = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)

        # Conversion from x,y,z to lon, lat, 0
        coords[:,0] = np.degrees(np.arctan2(coords[:,1], coords[:,0])) # longitude
        coords[:,1] = np.degrees(np.arcsin(coords[:,2] / R))           # latitude
        coords[:,2] = 0.0

        points.SetData(numpy_to_vtk(coords))

        # Get cell array
        cells         = ugrid.GetCells()
        cellformation = vtk_to_numpy(cells.GetData())
        num_cells     = ugrid.GetNumberOfCells()

        # Create point to cell mapping and an array of start indices of cells
        point2cells = [ [] for _ in range(num_points) ] # a list of lists
        startIDs    = [0] * num_cells # start index in cell data "cellformation"

        startID = 0
        for cell in range(num_cells) :              # loop through cells
            startIDs[cell] = startID
            size           = cellformation[startID] # number of vertices
            for i in range(size) :
                pid = cellformation[startID+1+i]
                point2cells[pid].append(cell)
                
            startID = startID + 1 + size

        # Loop over cells, check if it intersects with -180/180
        # separation line. If yes, move its vertices on one side horizontally
        # so that their lon coordinates are -180. Which side's vertices to be moved 
        # depends on numbers of vertices of that cell on both sides.
        if straight_edges : # straight edges
            for cell in range(num_cells) :    # loop through cells
                startID = startIDs[cell]
                size = cellformation[startID] # number of vertices

                # Fix with gaps near poles by shifting points near the poles to the poles
                lats = [i for i in range(size)]
                for i in range(size) :
                    lats[i]   = abs(coords[cellformation[startID+1+i],1])
                max_index = lats.index(max(lats))
                pid_max   = cellformation[startID+1+max_index]
                max_lat   = lats[max_index]
                if (max_lat > 80.0) :
                    coords[pid_max,1] =  np.sign(coords[pid_max,1]) * 90.0

                pids_pos = []
                pids_neg = []
                for i in range(size) :
                    pid = cellformation[startID+1+i]
                    if (  coords[pid,0] < -90.0 and coords[pid,0] > -180.0) : # vertex on negative side
                        pids_neg.append(pid)
                    elif (coords[pid,0] >  90.0 and coords[pid,0] <  180.0) : # vertex on positive side
                        pids_pos.append(pid)

                if (len(pids_pos) != 0 and len(pids_neg) != 0) : # cell intersects with -180/180 line
                    if (len(pids_pos) > len(pids_neg)) :         # move points on negative side
                        for pid in pids_neg :
                            coords[pid,0] = -180.0
                    else :                                       # move points on positive side
                        for pid in pids_pos :
                            coords[pid,0] = -180.0

            # Split points on the separation line
            last_pid = num_points
            for pid in range(num_points) :
                if (coords[pid,0] == -180.0 or coords[pid,0] == 180.0) :
                    coords[pid,0] = -180.0 
                    points.InsertNextPoint(180.0, coords[pid,1], coords[pid,2])
                    coords = np.vstack((coords, [[180.0, coords[pid,1], coords[pid,2]]]))

                    for cid in point2cells[pid] : # loop over all cells sharing point pid
                        # check cell is on positive or negative side
                        startID = startIDs[cid]
                        size = cellformation[startID] # number of vertices
                        is_pos = True
                        for i in range(size) :
                            p = cellformation[startID+1+i]
                            if (coords[p,0] != -180.0) :
                                if (coords[p,0]>0) :
                                    is_pos = True
                                    cellformation[startID+1+i] = last_pid
                                else :
                                    is_pos = False
                                break

                        # Replace pid with last_pid for those cells on positive-lon side
                        if (is_pos) :
                            loc = np.where(cellformation[startID+1:startID+1+size] == pid)
                            cells.ReplaceCellPointAtId(cid, int(loc[0][0]), last_pid)

                            startID = startIDs[cid]
                            size = cellformation[startID] # number of vertices
                            for i in range(size) :
                                if (cellformation[startID+1+i] == pid) :
                                    cellformation[startID+1+i] = last_pid
                                    break
                    last_pid += 1
                    
        else : # zig-zag
            points_on_sep = set()
            for cell in range(num_cells) :    # loop through cells
                startID  = startIDs[cell]
                size     = cellformation[startID] # number of vertices                
                pids_pos = []
                pids_neg = []
                for i in range(size) :
                    pid = cellformation[startID+1+i]
                    if (  coords[pid,0] < -365/5) : # on negative side
                        pids_neg.append(pid)
                    elif (coords[pid,0] >  365/5) : # on positive side
                        pids_pos.append(pid)

                if (len(pids_pos) > len(pids_neg)) : # positive triangle
                    points_on_sep.update(pids_neg)
                else :                               # negative triangle
                    points_on_sep.update(pids_pos)
                    startIDs[cell] = -startIDs[cell] # sign used to indicate which side cell should be placed on

            # Split points on separation line
            last_pid = num_points
            for pid in points_on_sep :
                # Add a new point
                if (coords[pid,0] > 0) :
                    new_lon = -360.0 + coords[pid,0]
                else :
                    new_lon =  360.0 + coords[pid,0]

                points.InsertNextPoint(new_lon, coords[pid,1], coords[pid,2])
                coords = np.vstack((coords, [[new_lon, coords[pid,1], coords[pid,2]]]))

                # Let newly added point be vertices of some triangles
                for cid in point2cells[pid] : # loop over all cells sharing point pid
                    # replace pid with last_pid for those cells sharing pid
                    startID = abs(startIDs[cid])
                    size = cellformation[startID]
                    if (startIDs[cid]>0 and coords[pid,0] < 0) or (startIDs[cid]<0 and coords[pid,0] > 0) : 
                        loc = np.where(cellformation[startID+1:startID+1+size] == pid)
                        loc  = int(loc[0][0])
                        cells.ReplaceCellPointAtId(cid, loc, last_pid)
                        cellformation[startID+1+loc] = last_pid
                last_pid += 1

        
        # Ensure no vertices are outside [-180,180] and that there are no gaps near -90, 90
        num_points = coords.shape[0]
        for pid in range(num_points) :
            if   (coords[pid,0] >  180.0) :
                coords[pid,0] = 180.0
            elif (coords[pid,0] < -180.0) :
                coords[pid,0] = -180.0

            if   (coords[pid,1] >  90.0 - 1.1*dtheta_max) :
                coords[pid,1] = 90.0
            elif (coords[pid,1] < -90.0 + 1.1*dtheta_max) :
                coords[pid,1] = -90.0

        # Update point data
        points.SetData(numpy_to_vtk(coords))

        if Delaunay=='y' :
            # Convert to point data
            cellToPoint = vtk.vtkCellDataToPointData()
            cellToPoint.SetInputData(ugrid)
            cellToPoint.Update()
            output = cellToPoint.GetOutput()

            # Apply vtkDelaunay2D filter
            delaunay = vtk.vtkDelaunay2D()
            delaunay.SetInputData(output)
            delaunay.Update()
            output = delaunay.GetOutput()

            point2cell = vtk.vtkPointDataToCellData()
            point2cell.SetInputData(output)
            point2cell.Update()
            output = point2cell.GetOutput()

            # Write out structured data
            writer = vtk.vtkXMLPolyDataWriter() 
            writer.SetFileName(outfile+".vtp")
            writer.SetInputData(output)
            writer.Write()
        else :
            writer = vtk.vtkUnstructuredGridWriter()
            writer.SetFileTypeToBinary()
            writer.SetFileName(outfile+".vtk")
            writer.SetInputData(ugrid)
            writer.Write()
            
for file in glob.glob("*tri_???_*.vtk"):
    os.remove(file)
