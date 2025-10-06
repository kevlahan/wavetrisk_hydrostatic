# Computes area integrated rms of either specified scalar or Rossby number in a zonal band  from longitude-latitude vtp files
import sys
import vtk
import glob
from utilities import *

rms_type = "scalar" # Rossby or scalar

Omega_planet = 7.29211e-5/6  # planet rotation

def rms_scalar (data, field, lat1, lat2) :
    # Integrated zonal rms statistics for field (e.g. Density) between latitudes lat1 and lat2

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
    scalar = band_with_area.GetCellData().GetArray(field)

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

        f = scalar.GetTuple1(i)
        num += (f*f) * A
        den += A

    # Area-weighted RMS
    return float(math.sqrt(num / den))

def rms_Rossby (data, lat1, lat2) :
    # Integrated zonal rms statistics for field (e.g. Density) between latitudes lat1 and lat2

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
    scalar   = band_with_area.GetCellData().GetArray("Vorticity")

    ncell = band_with_area.GetNumberOfCells()

    # --- cell centers to get latitude per cell (Y coord, degrees -> radians) ---
    centers = vtk.vtkCellCenters()
    centers.SetInputData(band_with_area)
    centers.Update()
    ctrs = centers.GetOutput()

    # --- accumulate area-weighted mean of Rossby^2 ---
    num = 0.0
    den = 0.0
    sin_eps = 1e-8
    for i in range(ncell):
        A = area_arr.GetTuple1(i)
        if A <= 0.0:
            continue

        # latitude (deg -> rad)
        lat_deg = ctrs.GetPoint(i)[1]
        phi = math.radians(lat_deg)
        s = math.sin(phi)
        if abs(s) < sin_eps:
            continue  # skip cells too close to equator (singular denominator)

        omega = scalar.GetTuple1(i)

        if not np.isfinite(omega):
            continue

        Ro = omega / (2.0 * Omega_planet * s)
        num += (Ro * Ro) * A
        den += A

    # Area-weighted RMS
    return float(math.sqrt(num / den))


#################################################################################################
#
#  Computes rms statistics of specified scalar or of Rossby number
#
#################################################################################################

# Input
if (len(sys.argv)<8) :
    print("\nUsage: python rms.py base_vtk_file k1 k2 t1 t2 dt field\n")
    print("Example: python3 rms.py drakeJ8Z60 1 60 120 120 5 Vorticity\n")
    print("run   = run prefix of vtp files to load (e.g. drakeJ8Z60)")
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
lat1  = float(sys.argv[7])
lat2  = float(sys.argv[8])
field = sys.argv[9] 

if rms_type == "Rossby":
    outfile = run+"_Rossby_rms.txt"
else:
    outfile = run+"_"+field+"_rms.txt"
    
f = open (outfile, "w")

for t in range (t1, t2+1):
    for k in range (k1, k2+1):
        infile = run+"_tri_lonlat_"+str(k).zfill(3)+"_"+str(t).zfill(4)+".vtp"
        print("rms %s of %s is" % (field, infile), end=" ")

        # Load the input vtk file
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(infile)
        reader.Update()
        data = reader.GetOutput()

        # Compute rms
        if rms_type == "Rossby":
            rms = rms_Rossby (data, lat1, lat2)
        elif rms_type == "scalar":
            rms = rms_scalar (data, field, lat1, lat2)
        else:
            print("Unknown rms type")
    
        print("%14.8e" % (rms)) 
        f.write("%10.4e %3d %14.8e\n" % (t*dt, k, rms))

    for file in glob.glob(run+'_tri_[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].vtk'):
        os.remove(file)



    

