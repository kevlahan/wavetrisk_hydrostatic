# Computes area integrated rms of specified vtk cell data
#
import sys
import vtk
from utilities import *

# Input
if (len(sys.argv)<5) :
    print("\nUsage: python rms.py base_vtk_file k1 k2 t1 t2 field\n")
    print("Example: python3 rms.py drakeJ8Z60_hex 54 54 120 120 5 density\n")
    print("base_vtk_file = prefix of vtk files to load (e.g. drakeJ8Z60_hex)")
    print("k1            = first vertical layer")
    print("k2            = last  vertical layer")
    print("t1            = first time count")
    print("t2            = last  time count")
    print("dt            = save interval (time = dt*count days)")
    print("field         = field to analyze (vorticity, density, OMEGA, etc\n")
    print("Output is saved to base_vtk_file.txt")
    exit(0)

file_base = sys.argv[1]
k1        = int(sys.argv[2])
k2        = int(sys.argv[3])
t1        = int(sys.argv[4])
t2        = int(sys.argv[5])
dt        = int(sys.argv[6])
field     = sys.argv[7]

outfile = file_base+"_rms.txt"
f = open(outfile, "w")

print()
for t in range (t1, t2+1):
    for k in range (k1, k2+1):
        infile = file_base+"_"+str(k).zfill(3)+"_"+str(t).zfill(4)+".vtk"
        print("rms %s of %s is" % (field, infile), end=" ")

        # Load the input vtk file
        vtk_data = vtk.vtkUnstructuredGridReader()
        vtk_data.ReadAllScalarsOn()
        vtk_data.SetFileName(infile)
        vtk_data.Update()

        # Compute rms
        rms = rms_int (vtk_data.GetOutput(), field)
    
        print("%14.8e" % (rms)) 
        f.write("%10.4e %3d %14.8e\n" % (t*dt, k, rms)) 

    


    

