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
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

verbose = False # echo file being processed

# Main program
if (len(sys.argv)<7):
    print("\nUsage: python xyz2lonlat.py run Jmin Jmax z1 z2 t1 t2 Delaunay\n")
    print("run      = file base name (without tri)")
    print("Jmin     = minimum level")
    print("Jmax     = maximum level")
    print("z1       = first z layer")
    print("z2       = last  z layer")
    print("t1       = first time")
    print("t2       = last time")
    print("Delaunay = y/n (interpolate to Delaunay grid to remove gaps and have straight edges)\n")
    print("output is file_lonlat run_tri_zzz_tttt.vtk or or run_tri_zzz_tttt.vtp (Delaunay = y)\n")
    print("Example:")
    print("python3 xyz2lonlat.py SimpleJ5J1Z30 5 7 1 30 0 100 y")
    exit(0)

run       = sys.argv[1]
Jmin      = int (sys.argv[2])
Jmax      = int (sys.argv[3])
z1        = int (sys.argv[4])
z2        = int (sys.argv[5])
t1        = int (sys.argv[6])
t2        = int (sys.argv[7])
Delaunay  = sys.argv[8]

# Grid dimensions (same number of rectangular cells as lozenge cells on the sphere) 
lat_dim = int (np.sqrt((10*4**Jmax + 2)/2))
lon_dim = 2*lat_dim

dtheta_min = 180/lat_dim
dtheta_max = dtheta_min * 2**(Jmax - Jmin)

for t in range (t1, t2+1):
    file = run+'_tri_'+str(t).zfill(4)+".vtk.tgz"
    if verbose: print ("Uncompressing file ", file)
    untar_files (file)
    file_vtk = run+"_tri_"+str(0).zfill(3)+"_"+str(t).zfill(4)+".vtk"
    if os.path.exists(file_vtk):
        os.remove(file_vtk)
    for z in range (z1, z2+1):
        # Load the input vtk file
        infile  = run+"_tri_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        outfile = run+"_tri_lonlat_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        if verbose: print ("Processing file ", infile)
        
        vtkreader = vtk.vtkDataSetReader()
        vtkreader.ReadAllScalarsOn()
        vtkreader.SetFileName(infile+".vtk")
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
                    if (np.abs(np.abs(coord[1])-90)<dtheta): # modify to avoid gaps if necessary
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
        writer.SetFileName(outfile+".vtp")
        writer.SetInputData(polydata)
        writer.Write()

    # Remove vtk files for time t
    pattern = f"{run}_tri_*_{str(t).zfill(4)}.vtk"
    for file in glob.glob(pattern):
        os.remove(file)
