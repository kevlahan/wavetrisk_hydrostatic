# Computes either averages or rms for given layers at given times. Returns mean and standard deviation over all times.
# of either specified scalar or Rossby number in a zonal band from longitude-latitude vtp files.
#
# Usage: python statistics.py stat_type base_vtk_file k1 k2 t1 t2 dt lat1 lat2 scl_Omega scl_radius field_type field
#
# field_type = Scalar, Rossby, KE, DeltaSM, DeltaI, VertFluxKE
#
# field (optional) selects field to analyze if field_type = Scalar
#
# Options for field:
#                    Level 
#                    Topography
#                    Penalization
#                    Ps/Eta
#                    Temperature/Density
#                    Velocity_Zonal
#                    Velocity_Meridional 
#                    OMEGA/Velocity_Vertical 
#                    Vorticity 
#                    Geopot_Height 
#                    P/Ps
#                    dz
#
# module load gcc mvapich python py-numpy py-scipy vtk
import sys
import vtk
import glob
import math
from utilities import *
from pathlib import Path

def compute_avg (data, field, lat1, lat2) :
    # Integrated zonal average statistics between latitudes lat1 and lat2

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

    if field_type == "Scalar":
        scalar = band_with_area.GetCellData().GetArray(field)
    
    zonal      = band_with_area.GetCellData().GetArray("Velocity_Zonal")
    meridional = band_with_area.GetCellData().GetArray("Velocity_Meridional")
    vertical   = band_with_area.GetCellData().GetArray("Velocity_Vertical")
    vorticity  = band_with_area.GetCellData().GetArray("Vorticity")

    ncell = band_with_area.GetNumberOfCells()

    # --- cell centers to get latitude per cell (Y coord, degrees -> radians) ---
    centers = vtk.vtkCellCenters()
    centers.SetInputData(band_with_area)
    centers.Update()
    ctrs = centers.GetOutput()

    # --- accumulate area-weighted mean ---
    num = 0.0
    den = 0.0
    for i in range(ncell):
        A = area_arr.GetTuple1(i)
        if A <= 0.0:
            continue
        den += A
        
        theta = math.radians(ctrs.GetPoint(i)[0]) # longitude 
        phi   = math.radians(ctrs.GetPoint(i)[1]) # latitude

        if field_type == "Scalar":
            sclr = scalar.GetTuple1(i)
            
            num += sclr * A
        elif field_type == "Rossby":
            omega = vorticity.GetTuple1(i)
            f     = 2.0 * Omega_planet * math.sin(phi)
            Ro    = abs(omega / f)
            
            num += Ro * A
        elif field_type == "KE":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            
            num += 0.5 * math.sqrt(u * u + v * v) * A
        elif field_type == "DeltaSM":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            f     = 2.0 * Omega_planet * math.sin(phi)
            
            num +=  math.sqrt(u * u + v * v) / f * A
        elif field_type == "DeltaI":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            beta  = 2.0 * Omega_planet * math.cos(phi) / Radius_planet
            
            num +=  math.sqrt(math.sqrt(u * u + v * v) / beta) * A
        elif field_type == "VertFluxKE":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            w     = vertical.GetTuple1(i)
            KE    = 0.5 * (u * u + v * v)
            
            num += KE * w * A
    
    # Area-weighted average
    return num / den

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
    
    if field_type == "Scalar":
        scalar = band_with_area.GetCellData().GetArray(field)
        
    zonal      = band_with_area.GetCellData().GetArray("Velocity_Zonal")
    meridional = band_with_area.GetCellData().GetArray("Velocity_Meridional")
    vertical   = band_with_area.GetCellData().GetArray("Velocity_Vertical")
    vorticity  = band_with_area.GetCellData().GetArray("Vorticity")

    ncell = band_with_area.GetNumberOfCells()

    # --- cell centers to get latitude per cell (Y coord, degrees -> radians) ---
    centers = vtk.vtkCellCenters()
    centers.SetInputData(band_with_area)
    centers.Update()
    ctrs = centers.GetOutput()

    # --- accumulate area-weighted means ---
    num = 0.0
    den = 0.0
    for i in range(ncell):
        A = area_arr.GetTuple1(i)
        if A <= 0.0:
            continue
        den += A
        
        theta = math.radians(ctrs.GetPoint(i)[0]) # longitude 
        phi   = math.radians(ctrs.GetPoint(i)[1]) # latitude

        if field_type == "Scalar":
            sclr = scalar.GetTuple1(i)
            
            num += (sclr * sclr) * A
        elif field_type == "Rossby":
            omega = vorticity.GetTuple1(i)
            f     = 2.0 * Omega_planet * math.sin(phi)
            Ro    = omega / f
            
            num += (Ro * Ro) * A
        elif stat_type == "KE":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            KE    = 0.5 * (u * u + v * v)

            num += (KE * KE) * A
        elif stat_type == "DeltaSM":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            f     = 2.0 * Omega_planet * math.sin(phi)
            
            num +=  (u * u + v * v) / (f * f) * A
        elif stat_type == "DeltaI":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            beta  = 2.0 * Omega_planet * math.cos(phi) / Radius_planet
            
            num +=  math.sqrt(u * u + v * v) / beta * A
        elif avg_type == "VertFluxKE":
            u     = zonal.GetTuple1(i)
            v     = meridional.GetTuple1(i)
            w     = vertical.GetTuple1(i)
            KE    = 0.5 * (u * u + v * v)
            
            num += (KE * w) * (KE * w)  * A
    
    # Area-weighted RMS
    return float(math.sqrt(num / den))

#################################################################################################
#
#  Computes averages statistics of specified scalar or other statistical quantities
#
#################################################################################################

# Input
if (len(sys.argv)<10) :
    print("Computes either averages or rms for given layers at given times. Returns mean and standard deviation over all times.")
    print("\nUsage: python statistics.py avg_type base_vtk_file k1 k2 t1 t2 dt lat1 lat2 scl_Omega scl_radius field_type field\n")
    print("Example 1: python statistics.py avg drakeJ8Z60 1 60 120 120 5 15 50 6 6 Rossby \n")
    print("Example 2: python statistics.py rms drakeJ8Z60 1 60 120 120 5 15 50 6 6 Scalar Vorticity \n")
    print("stat_type  = avg or rms")
    print("run        = run prefix of vtp files to load (e.g. drakeJ8Z60)")
    print("k1         = First vertical layer")
    print("k2         = Last  vertical layer")
    print("t1         = First time count")
    print("t2         = Last  time count")
    print("dt         = Save interval (time = dt*count days)")
    print("lat1       = Minimum latitude")
    print("lat2       = Maximum latitude")
    print("scl_Omega  = scaling factor for Omega:  Omega_Earth /scl_Omega")
    print("scl_radius = scaling factor for radius: radius_Earth/scl_radius")
    print("field_type = Scalar, Rossby, speed, deltaSM, deltaI")
    print("field      = Field to analyze: \n \
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
    print("Output is saved to run_[field_type].txt")
    exit(0)

stat_type  = sys.argv[1]
run        = sys.argv[2]
k1         = int(sys.argv[3])
k2         = int(sys.argv[4])
t1         = int(sys.argv[5])
t2         = int(sys.argv[6])
dt         = int(sys.argv[7])
lat1       = float(sys.argv[8])
lat2       = float(sys.argv[9])
scl_Omega  = float(sys.argv[10])
scl_radius = float(sys.argv[11])
field_type = sys.argv[12]

if field_type == "Scalar":
    field = sys.argv[13] 

if field_type == "Scalar":
    outfile = run+"_"+field+"_"+stat_type+".txt"
else:
    field = ""
    outfile = run+'_'+field_type+"_"+stat_type+".txt"

Omega_planet  = 7.29211e-5 / scl_Omega  # planet rotation
Radius_planet = 6371.229e3 / scl_radius # planet radius
    
f = open (outfile, "w")

import math

for k in range(k1, k2+1):
    ntime = 0
    mean  = 0.0      # running mean
    M2    = 0.0      # running sum of squares of differences from the mean

    for t in range(t1, t2+1):
        infile = run+"_tri_lonlat_"+str(k).zfill(3)+"_"+str(t).zfill(4)+".vtp"
        print(f"Processing file: {infile}", flush=True)

        p = Path(infile).expanduser()
        if not p.is_file():
            sys.exit(f"ERROR: file not found: {p}")   # exits with code 1

        # Load the input vtk file
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(infile)
        reader.Update()
        data = reader.GetOutput()

        # Value at this time
        if stat_type == "avg":
            x = compute_avg(data, field, lat1, lat2)
        elif stat_type == "rms":
            x = compute_rms(data, field, lat1, lat2)

        # --- Welford update ---
        ntime += 1
        delta  = x - mean
        mean  += delta / ntime
        delta2 = x - mean
        M2    += delta * delta2
        # ----------------------

    # After the loop over t, compute variance / std

    if ntime > 1:
        variance = M2 / (ntime - 1) # sample variance (divide by n-1)
    else:
        variance = 0.0

    # Guard against tiny negative due to rounding
    std_t = math.sqrt(max(variance, 0.0))

    avg_t = mean

    # Save: k, mean, std deviation
    f.write("%3d %14.8e %14.8e\n" % (k, avg_t, std_t))

# Clean up temporary .vtk files (once, at the end)
for file in glob.glob(run+'_tri_[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].vtk'):
    os.remove(file)


