import numpy as np
import sys
import vtk

# Area/volume integrated cell data example

# Input
if (len(sys.argv)<1) :
    print("Usage: python area.py vtk_file")
    exit(0)

infile = sys.argv[1]

# Load the input vtk file
vtk_data = vtk.vtkUnstructuredGridReader()
vtk_data.ReadAllScalarsOn()
vtk_data.SetFileName(infile)
vtk_data.Update()

integrate = vtk.vtkIntegrateAttributes()
integrate.SetInputData(vtk_data.GetOutput())
integrate.SetDivideAllCellDataByVolume(1)
integrate.Update()

integrated = integrate.GetOutput()

cd = integrated.GetCellData()

print("Integrated cell data")
for i in range(cd.GetNumberOfArrays()):
    arrayName = cd.GetArray(i).GetName()
    arrayValue = cd.GetArray(i).GetTuple(0)[0]
    print("{} = {:e}".format(arrayName, arrayValue))
