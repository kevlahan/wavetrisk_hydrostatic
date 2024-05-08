# Converts xyz data in a vtk file to lonlat data and saves in a new file
# (provides a simple longitude-latitude projection of spherical data)
# Usage: python xyz2lonlat.py input_vtk_file output_vtk_file
# Weiguang Guan (SHARCNET) 2024-04-16

import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# Split a cell specifyed by "cellid" in an unstructured grid specified by "ugrid"
def split_cell(cellid, offset, ugrid) :
    points = ugrid.GetPoints()
    coords = vtk_to_numpy(points.GetData())
    cells = ugrid.GetCells().GetData()
    cellformation = vtk_to_numpy(cells)
    celldata = ugrid.GetCellData()
    numpoints = points.GetNumberOfPoints()

    size = cellformation[offset]
    for i in range(size) : # loop through vertices of the cell to be split
        pid = cellformation[offset+1+i]
        #print(coords[pid,:])
        if (coords[pid,0]>0) :
            points.InsertNextPoint(coords[pid,:])
            points.SetPoint(pid, -180.0, coords[pid,1], coords[pid,2])
        else :
            points.InsertNextPoint(180.0, coords[pid,1], coords[pid,2])

    # Construct/insert a cell with newly inserted points
    point_ids = list(range(numpoints, (numpoints+size)))
    ugrid.InsertNextCell(vtk.VTK_POLYGON , len(point_ids), point_ids)

    numattributes = celldata.GetNumberOfArrays()
    for i in range(numattributes) :
        attrarray = celldata.GetArray(i)
        attrarray.InsertNextTuple1(attrarray.GetTuple1(cellid))


# Main program
if (len(sys.argv)<3) :
    print("Usage: python xyz2lonlat.py input_vtk_file output_vtk_file")
    exit(0)

infile = sys.argv[1]
outfile = sys.argv[2]

# Load the input vtk file
vtkreader = vtk.vtkUnstructuredGridReader()
vtkreader.ReadAllScalarsOn()
vtkreader.SetFileName(infile)
vtkreader.Update()

# Get the unstructed grid data
unstrctGrid = vtkreader.GetOutput()

# Get coordinates of vertices
points = unstrctGrid.GetPoints()
coords = vtk_to_numpy(points.GetData())

# Conversion from x,y,z to log, lat, 0
R = 6.3707e6

coords[:,0] = np.degrees(np.arctan2(coords[:,1], coords[:,0])) # longitude
coords[:,1] = np.degrees(np.arcsin(coords[:,2] / R)) # latitude
coords[:,2] = 0.0

points.SetData(numpy_to_vtk(coords))

# Re-adjust coords
cellformation = np.copy(vtk_to_numpy(unstrctGrid.GetCells().GetData()))

num_cells = unstrctGrid.GetNumberOfCells()
startID = 0
for cell in range(num_cells) : # loop through cells
    size = cellformation[startID] #cells.GetCellSize(cell)  # number of vertices

    # Check if it is a cell on the edge
    num_positives = 0
    for i in range(size) : # loop through vertices of a cell
        pid = cellformation[startID+1+i]
        if (coords[pid,0]>0) :
            num_positives += 1

    if not (num_positives==size or num_positives==0) : # perhaps needs to split
        if (abs(coords[cellformation[startID+1],0])>90) :
            split_cell(cell, startID, unstrctGrid)

    startID = startID + 1 + size

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileTypeToBinary()
writer.SetFileName(outfile)
writer.SetInputData(unstrctGrid)
writer.Write()
