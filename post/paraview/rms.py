# Computes rms of vtk cell data
#
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# Define rms function
def rms(arr):
    return np.sqrt(np.mean(np.square(arr)))

# Main program
if (len(sys.argv)<5) :
    print("Usage: python energy.py base_vtk_file tstart tend field")
    print("base_vtk_file = prefix of vtk files to load (e.g. HS_J7_001)")
    print("tstart = first file")
    print("tend   = last last file")
    print("dt     = save interval (days)")
    print("field  = field to analyze (mass, vorticity, etc)")
    print("Output saved to base_vtk_file.txt")
    exit(0)

file_base = sys.argv[1]
j1        = int(sys.argv[2])
j2        = int(sys.argv[3])
dt        = int(sys.argv[4])
field     = sys.argv[5]

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
    data = cell_data.GetArray(field)

    rms_data = rms(data)

    print("{:>2d} {:e}".format(dt*j, rms_data))
    f.write("{:>2d} {:e}\n".format(dt*j, rms_data))

    # num_cells = vtk_data.GetOutput().GetNumberOfCells()
    
    # calc = vtk.vtkArrayCalculator()
    # calc.SetInputData(vtk_data.GetOutput())
    # calc.SetFunction("mass**2")
    # calc.SetResultArrayName("square")
    # calc.Update()

    # integrate = vtk.vtkIntegrateAttributes()
    # integrate.SetInputData(vtk_data.GetOutput())
    # integrate.SetDivideAllCellDataByVolume(1)
    # integrate.Update()

    # integrated = integrate.GetOutput()
    # cd = integrated.GetCellData()

    # print("Integrated cell data")
    # for i in range(cd.GetNumberOfArrays()) :
    #     arrayName = cd.GetArray(i).GetName()
    #     arrayValue = cd.GetArray(i).GetTuple(0)[0]
    #     print("%s: %f" % (arrayName, arrayValue))


    

