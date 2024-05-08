# Computes rms of vtk cell data
#
import sys
import numpy as np
import vtk



# Define rms function
def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

# Main program
if (len(sys.argv)<4) :
    print("Usage: python energy.py base_vtk_file tstart tend field")
    print("tstart = first file")
    print("tend   = last last file")
    print("field  = field to analyze (mass, vorticity, etc)")
    exit(0)

file_base = sys.argv[1]
j1        = int(sys.argv[2])
j2        = int(sys.argv[3])
field     = sys.argv[4]

outfile = file_base+"_rms.txt"
f = open(outfile, "w")

for j in range (j1, j2+1):
    infile = file_base+".1"+str(j).zfill(4)+".vtk"

    # Load the input vtk file
    vtk_data = vtk.vtkUnstructuredGridReader()
    vtk_data.ReadAllScalarsOn()
    vtk_data.SetFileName(infile)
    vtk_data.Update()

    cell_data = vtk_data.GetOutput().GetCellData()

    data      = cell_data.GetArray(field)

    calc = vtk.vtkArrayCalculator()
    calc.SetInputData(vtk_data.GetOutput())
    calc.SetFunction("5")
    calc.SetResultArrayName("MyResults")
    calc.Update()

    num_cells = vtk_data.GetOutput().GetNumberOfCells()
    print(vtk)

    #calculator = vtk.vtkCalculator()

    integrate = vtk.vtkIntegrateAttributes()

    integrate.SetInputData(vtk_data.GetOutput())
    integrate.SetDivideAllCellDataByVolume(1)
    integrate.Update()

    integrated = integrate.GetOutput()
    cd = integrated.GetCellData()

    print("Integrated cell data")
    for i in range(cd.GetNumberOfArrays()) :
        arrayName = cd.GetArray(i).GetName()
        arrayValue = cd.GetArray(i).GetTuple(0)[0]
        print("%s: %f" % (arrayName, arrayValue))

    
    f.write("{:d} {:e}\n".format(j, rms(data)))

