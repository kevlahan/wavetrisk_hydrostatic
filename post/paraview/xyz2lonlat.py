import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

if (len(sys.argv)<3) :
    print("Usage: python nick.py input_vtk_file output_vtk_file")
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
print(coords.shape)

# Conversion from x,y,z to log, lat, 0
R = 6.3707e6

coords[:,0] = np.degrees(np.arctan2(coords[:,1], coords[:,0])) # longitude
coords[:,1] = np.degrees(np.arcsin(coords[:,2] / R)) # latitude
coords[:,2] = 0.0

# Re-adjust coords
cells = unstrctGrid.GetCells()
cellformation = vtk_to_numpy(cells.GetData())

num_cells = unstrctGrid.GetNumberOfCells()
startID = 0
for cell in range(num_cells) : # loop through cells
    #print("Cell:", cell)
    size = cellformation[startID] #cells.GetCellSize(cell)  # number of vertices

    # Check if it is a cell on the edge
    num_positives = 0
    for pid in range(size) : # loop through vertices of a cell
        if (coords[cellformation[startID+1+pid],0]>0) :
            num_positives += 1

    if not (num_positives==size or num_positives==0) : # need to adjust
        #print("Cell to adjust:", cell)
        if (num_positives>size/2) : # Make all positive
            for pid in range(size) :
                if (coords[cellformation[startID+1+pid],0]<-90) :
                    coords[cellformation[startID+1+pid],0] += 360.0
        else : # Make all negative
            for pid in range(size) :
                if (coords[cellformation[startID+1+pid],0]>90) :
                    coords[cellformation[startID+1+pid],0] -= 360.0

    startID = startID + 1 + size

# Put new coords (lon-lat) back into the vtkUnstructuredGrid object
points.SetData(numpy_to_vtk(coords))
#unstrctGrid.SetPoints(points)

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileTypeToBinary()
writer.SetFileName(outfile)
writer.SetInputData(unstrctGrid)
writer.Write()
