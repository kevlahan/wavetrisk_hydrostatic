# Computes area integrated rms of either specified scalar or Rossby number in a zonal band  from longitude-latitude vtp files
# Usage: python rms.py base_vtk_file k1 k2 t1 t2 dt lat1 lat2 rms_type field
# rms_type = scalar, Rossby, speed, deltaSM, deltaI
# field    = Field to analyze:
#                Options =
#                   Level 
#                   Topography
#                   Penalization
#                   Ps/Eta
#                   Temperature/Density
#                   Velocity_Zonal
#                   Velocity_Meridional 
#                   OMEGA/Velocity_Vertical 
#                   Vorticity 
#                   Geopot_Height 
#                   P/Ps
#                   dz
import sys
import vtk
import glob
from utilities import *

Omega_planet  = 7.29211e-5 / 6  # planet rotation
Radius_planet = 6371.229e3 / 6  # planet radius

def compute_rms (data, field, lat1, lat2) :
    # Integrated zonal rms statistics between latitudes lat1 and lat2

    # --- normalize zonal region ---
    if lat1 > lat2:
        lat1, lat2 = lat2, lat1

    # --- extract zonal region ---
    xmin, xmax, ymin, ymax, zmin, zmax = data.GetBounds()
    lo, hi = max(lat1, ymin), min(lat2, ymax)

    box = vtk.vtkBox()
    box.SetBounds(xmin, xmax, lo, hi, zmin, zmax)

    ext = vtk.vtkExtractPolyDataGeometry()
    ext.SetInputData(data)
    ext.SetImplicitFunction(box)
    ext.ExtractInsideOn()
    ext.ExtractBoundaryCellsOn()
    ext.Update()
    band = ext.GetOutput()

    # --- compute per-cell areas (works for polygonal polydata) ---
    cellsz = vtk.vtkCellSizeFilter()
    cellsz.SetInputData(band)
    cellsz.SetComputeArea(True)
    cellsz.SetComputeVolume(False)
    cellsz.SetComputeLength(False)
    cellsz.Update()
    band_with_area = cellsz.GetOutput()

    area_arr = band_with_area.GetCellData().GetArray("Area")
    
    scalar     = band_with_area.GetCellData().GetArray(field)
    zonal      = band_with_area.GetCellData().GetArray("Velocity_Zonal")
    meridional = band_with_area.GetCellData().GetArray("Velocity_Meridional")
    vorticity  = band_with_area.GetCellData().GetArray("Vorticity")

    ncell = band_with_area.GetNumberOfCells()

    # --- cell centers to get latitude per cell (Y coord, degrees -> radians) ---
    centers = vtk.vtkCellCenters()
    centers.SetInputData(band_with_area)
    centers.Update()
    ctrs = centers.GetOutput()

    # --- accumulate area-weighted mean of Rossby^2 ---
    num = 0.0
    den = 0.0
    for i in range(ncell):
        A = area_arr.GetTuple1(i)
        if A <= 0.0:
            continue
        den += A
        
        theta = math.radians(ctrs.GetPoint(i)[0]) # longitude 
        phi   = math.radians(ctrs.GetPoint(i)[1]) # latitude

        if rms_type == "scalar":
            sclr = scalar.GetTuple1(i)
            
            num += (sclr * sclr) * A
        elif rms_type == "Rossby":
            omega = vorticity.GetTuple1(i)
            f     = (2.0 * Omega_planet * math.sin(phi))
            Ro    = omega / f
            
            num += (Ro * Ro) * A
        elif rms_type == "speed":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            speed = math.sqrt(u * u + v * v)
            
            num += speed * A
        elif rms_type == "deltaSM":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            f     = 2.0 * Omega_planet * math.sin(phi)
            speed = math.sqrt(u * u + v * v)
            
            num +=  (speed / f) * A
        elif rms_type == "deltaI":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            speed = math.sqrt(u * u + v * v)
            beta  = 2.0 * Omega_planet * math.cos(phi) / Radius_planet
            
            num += math.sqrt(speed / beta) * A
    
    # Area-weighted RMS
    return float(math.sqrt(num / den))


#################################################################################################
#
#  Computes rms statistics of specified scalar or of Rossby number
#
#################################################################################################

# Input
if (len(sys.argv)<10) :
    print("\nUsage: python rms.py base_vtk_file k1 k2 t1 t2 dt lat1 lat2 rms_type field\n")
    print("Example 1: python3 rms.py drakeJ8Z60 1 60 120 120 5 Rossby \n")
    print("Example 2: python3 rms.py drakeJ8Z60 1 60 120 120 5 scalar Vorticity \n")
    print("run      = run prefix of vtp files to load (e.g. drakeJ8Z60)")
    print("k1       = First vertical layer")
    print("k2       = Last  vertical layer")
    print("t1       = First time count")
    print("t2       = Last  time count")
    print("dt       = Save interval (time = dt*count days)")
    print("lat1     = Minimum latitude")
    print("lat2     = Maximum latitude")
    print("rms_type = scalar, Rossby, speed, deltaSM, deltaI")
    print("field    = Field to analyze: \n \
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

run      = sys.argv[1]
k1       = int(sys.argv[2])
k2       = int(sys.argv[3])
t1       = int(sys.argv[4])
t2       = int(sys.argv[5])
dt       = int(sys.argv[6])
lat1     = float(sys.argv[7])
lat2     = float(sys.argv[8])
rms_type = sys.argv[9]
if rms_type == "scalar":
    field = sys.argv[10] 

if rms_type == "scalar":
    outfile = run+"_"+field+"_rms.txt"
elif rms_type == "Rossby":
    field = ""
    outfile = run+"_Rossby_rms.txt"
elif rms_type == "speed":
    field = ""
    outfile = run+"_speed_rms.txt"
elif rms_type == "deltaSM":
    field = ""
    outfile = run+"_deltaSM_rms.txt"
elif rms_type == "deltaI":
    field = ""
    outfile = run+"_deltaI_rms.txt"
    
f = open (outfile, "w")

for t in range (t1, t2+1):
    for k in range (k1, k2+1):
        infile = run+"_tri_lonlat_"+str(k).zfill(3)+"_"+str(t).zfill(4)+".vtp"

        if rms_type == "scalar":
            print("rms %s of %s is" % (field, infile), end=" ")
        else:
            print("rms %s of %s is" % (rms_type, infile), end=" ")

        # Load the input vtk file
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(infile)
        reader.Update()
        data = reader.GetOutput()

        # Compute rms
        rms = compute_rms (data, field, lat1, lat2)

        print("%14.8e" % (rms)) 
        f.write("%10.4e %3d %14.8e\n" % (t*dt, k, rms))

    for file in glob.glob(run+'_tri_[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].vtk'):
        os.remove(file)



    

