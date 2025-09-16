# Computes area integrated rms of specified vtk cell data
import sys
import vtk
import glob
from utilities import *

# Input
if (len(sys.argv)<8) :
    print("\nUsage: python rms.py base_vtk_file k1 k2 t1 t2 dt field\n")
    print("Example: python3 rms.py drakeJ8Z60 1 60 120 120 5 Vorticity\n")
    print("run   = run prefix of vtk files to load (e.g. drakeJ8Z60)")
    print("k1    = First vertical layer")
    print("k2    = Last  vertical layer")
    print("t1    = First time count")
    print("t2    = Last  time count")
    print("dt    = Save interval (time = dt*count days)")
    print("field = Field to analyze: \n \
               Options = \n \
                  Level \n \
                  Topography \n \
                  Penalization \n \
                  Ps/Eta \n \
                  Temperature/Density \n \
                  Velocity_Zonal \n \
                  Velocity_Meridional \n \
                  OMEGA/Velocity_Vertical \n \
                  Vorticity \n \
                  Geopot_Height \n \
                  P/Ps \n \
                  dz \n")
    print("Output is saved to base_vtk_file.txt")
    exit(0)

run   = sys.argv[1]
k1    = int(sys.argv[2])
k2    = int(sys.argv[3])
t1    = int(sys.argv[4])
t2    = int(sys.argv[5])
dt    = int(sys.argv[6])
field = sys.argv[7] 

outfile = run+"_"+field+"_rms.txt"
f = open (outfile, "w")

for t in range (t1, t2+1):
    file = run+'_tri_'+str(t).zfill(4)+".vtk.tgz"
    untar_files (file)
    for k in range (k1, k2+1):
        infile = run+"_tri_"+str(k).zfill(3)+"_"+str(t).zfill(4)+".vtk"
        print("rms %s of %s is" % (field, infile), end=" ")

        # Load the input vtk file
        vtk_data = vtk.vtkDataSetReader()
        vtk_data.ReadAllScalarsOn()
        vtk_data.SetFileName(infile)
        vtk_data.Update()

        # Compute rms
        rms = rms_int (vtk_data.GetOutput(), field)
    
        print("%14.8e" % (rms)) 
        f.write("%10.4e %3d %14.8e\n" % (t*dt, k, rms)) 

    for file in glob.glob(run+'_tri_[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].vtk'):
        os.remove(file)


    

