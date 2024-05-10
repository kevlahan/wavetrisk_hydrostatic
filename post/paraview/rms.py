# Computes area integrated rms of specified vtk cell data
#
import sys
import vtk
from utilities import *

# Input
if (len(sys.argv)<5) :
    print("Usage: python energy.py base_vtk_file tstart tend dt field")
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

    # Compute rms
    rms = rms_int (vtk_data.GetOutput(), "mass")
    
    print("{:e}".format(rms))
    f.write("{:>4d} {:e}\n".format(dt*j, rms))

    


    

