# Computes either averages or rms for given layers at given times. Returns mean and standard deviation over all times.
# of either specified scalar or Rossby number in a zonal band from longitude-latitude vtp files.
#
# Example: python statistics.py stat_type base_vtk_file k1 k2 t1 t2 dt lat1 lat2 scl_Omega scl_radius field_type field
#
# field_type = Scalar, Rossby, KE, DeltaSM, DeltaI, VertFluxKE, EnergyDiagnostics
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
#                    Vertical_Velocity/Velocity_Vertical 
#                    Vorticity 
#                    Geopot_Height 
#                    P/Ps
#                    dz
#
#
# EnergyDiagnostics: produces five terms in one labelled twelve-column table.
#
# Example: python statistics.py avg drakeJ8Z60 1 60 120 120 5 15 50 6 6 EnergyDiagnostics
#
# Supply the complete water column, ordered bottom-to-top by default.
#
# Optional:
# --layer-order top-to-bottom --g 9.80616 --rho-0 1030
# --pressure-density in-situ (default: Density*g*dz); reference is a comparison option.
# --surface-pressure 0 sets surface gauge pressure in Pa.
# --height-field Geopot_Height uses exported physical centre z for derivatives;
# otherwise heights are reconstructed from dz below a flat surface z=0.

# Output organization: 12 columns for each layer (includes time mean and std deviation for each quantity)  
#
# Column	Quantity
# 1	        Layer number                              k
# 2      	Time-mean centre height                   z           (m)
# 3–4   	Eddy buoyancy conversion                <b'w'>        (W/kg)
# 5–6	        Eddy pressure-energy flux               <p'w'>        (W/m²)
# 7–8	        Eddy pressure-flux convergence     -d_z <p'w'> / rho0 (W/kg)
# 9–10	        Eddy KE vertical advective flux   rho_0 <w K>         (W/m²)
# 11–12	        Eddy KE advective-flux convergence -d_z <w K>         (W/kg)
#
# module load gcc mvapich python py-numpy py-scipy vtk 
import sys
import os
import argparse
import numpy as np
import vtk
import glob
import math
from utilities import *
from pathlib import Path

def get_vertical_velocity(cells):
    """Physical incompressible vertical velocity, never pressure-coordinate OMEGA."""
    vertical = cells.GetArray("Vertical_Velocity")
    return vertical if vertical is not None else cells.GetArray("Velocity_Vertical")


# Essential diagnostics. At each layer/time, <q> is the area-weighted
# latitude-band mean and q'=q-<q>; primes are spatial, not temporal.
# b=-g*(Density-rho_0)/rho_0; p_d=p+rho_0*g*z; z is positive upward.
#
# Output: layer, time-mean centre height (m), then mean/std pairs for each energy flux term.
#
# Time std is the sample standard deviation across snapshots (zero for one).
ENERGY_TERMS = [
    ("HeightCentre", "m", "<z>, positive upward"), 
    
    ("EddyBuoyancyConversion", "W/kg",
     "<w'b'>; potential-to-eddy-kinetic energy conversion; positive is kinetic-energy gain"),
    
    ("EddyPressureEnergyFlux", "W/m^2",
     "<p_d'w'>; vertical eddy pressure-energy transport; positive is upward"),
    
    ("EddyPressureFluxProfileConvergence", "W/kg",
     "-d_z(<p_d'w'>)/rho_0; kinetic-energy gain from pressure-flux convergence; positive is gain"),
    
    ("EddyKineticEnergyAdvectiveFlux", "W/m^2",
     "rho_0*<w*K_e>, K_e=(u'^2+v'^2)/2; vertical advection of eddy horizontal KE by full w; positive is upward"),
    
    ("EddyKineticEnergyFluxProfileConvergence", "W/kg",
     "-d_z(<w*K_e>); eddy KE gain from vertical advective-flux convergence; positive is gain"),
]


def reconstruct_hydrostatic(dz, density, gravity, reference_density,
                            pressure_density="in-situ", surface_pressure=0.0):
    """Top-to-bottom arrays. Integrate upward then impose p_surface.

    dp/dz=-rho*g. Each column has its own pressure integration constant.
    Density is piecewise constant within each layer. Returns centre and edge p,
    and centre z relative to a flat surface at z=0. No bottom p=0 assumption.
    """
    dz = np.asarray(dz, dtype=float)
    density = np.asarray(density, dtype=float)
    if dz.ndim != 2 or dz.shape != density.shape:
        raise ValueError("dz and Density must be matching layer-by-column arrays")
    if not np.all(np.isfinite(dz)) or np.any(dz <= 0):
        raise ValueError("EnergyDiagnostics requires finite positive dz in every selected column")
    if not np.all(np.isfinite(density)) or np.any(density <= 0):
        raise ValueError("EnergyDiagnostics requires finite positive Density")
    rho = density if pressure_density == "in-situ" else reference_density
    jumps = gravity * rho * dz
    edges = np.zeros((dz.shape[0] + 1, dz.shape[1]))
    z_edges = np.zeros_like(edges)
    for k in range(dz.shape[0] - 1, -1, -1):
        edges[k] = edges[k + 1] - jumps[k]
        z_edges[k] = z_edges[k + 1] + dz[k]
    edges += surface_pressure - edges[0].copy()
    z_edges -= z_edges[0].copy()
    return (edges[:-1] + edges[1:]) / 2, edges, (z_edges[:-1] + z_edges[1:]) / 2


def energy_profile(dz, density, u, v, w, area, options, heights=None):
    """Return mean height and the five instantaneous diagnostics.

    Pressure retains full density by default. Its profile convergence uses
    mean layer heights; it is not a conservative moving-layer budget term.
    """
    pressure, edges, z = reconstruct_hydrostatic(
        dz, density, options.g, options.rho_0,
        options.pressure_density, options.surface_pressure)
    if heights is not None:
        z = heights
    weights = area / np.sum(area, axis=1, keepdims=True)
    def average(x):
        return np.sum(weights * x, axis=1)
    def anomaly(x):
        return x - average(x)[:, None]
    buoyancy = -options.g * (density - options.rho_0) / options.rho_0
    eddy_b = average(anomaly(w) * anomaly(buoyancy))
    reduced_pressure = pressure + options.rho_0 * options.g * z
    eddy_p = average(anomaly(w) * anomaly(reduced_pressure))
    mean_z = average(z)
    convergence = -vertical_profile_derivative(eddy_p, mean_z) / options.rho_0
    # Eddy horizontal KE is advected by the full vertical velocity w.
    # This includes both mean-flow and fluctuating vertical advection.
    eddy_ke = 0.5 * (anomaly(u)**2 + anomaly(v)**2)
    ke_flux = options.rho_0 * average(w * eddy_ke)
    ke_convergence = -vertical_profile_derivative(ke_flux, mean_z) / options.rho_0
    return np.column_stack([mean_z, eddy_b, eddy_p, convergence,
                            ke_flux, ke_convergence])


def write_energy_diagnostics(output, run, k1, k2, t1, t2, lat1, lat2, options):
    """Read a complete surface-to-bottom column; preserve original layer labels."""
    layers = list(range(k1, k2 + 1))
    if options.layer_order == "bottom-to-top":
        layers.reverse()
    mean = np.zeros((len(layers), len(ENERGY_TERMS)))
    M2 = np.zeros_like(mean)
    for count, t in enumerate(range(t1, t2 + 1), 1):
        fields = {name: [] for name in ("dz", "Density", "u", "v", "w", "Area", "z")}
        reference_xy = None
        for k in layers:
            path = Path(f"{run}_tri_lonlat_{k:03d}_{t:04d}.vtp").expanduser()
            print(f"Processing file: {path}", flush=True)
            if not path.is_file():
                raise FileNotFoundError(f"file not found: {path}")
            reader = vtk.vtkXMLPolyDataReader(); reader.SetFileName(str(path)); reader.Update()
            band = extract_band_with_area(reader.GetOutput(), lat1, lat2)
            cells = band.GetCellData()
            def values(name):
                a = get_vertical_velocity(cells) if name == "w" else cells.GetArray(name)
                if a is None or a.GetNumberOfComponents() != 1:
                    raise ValueError(f"{path}: missing scalar {name} (w requires Vertical_Velocity)")
                vals = np.array([a.GetTuple1(i) for i in range(band.GetNumberOfCells())])
                if not np.all(np.isfinite(vals)):
                    raise ValueError(f"{path}: non-finite {name}")
                return vals
            area = values("Area")
            active = area > 0
            if not np.any(active):
                raise ValueError(f"{path}: no positive-area cells in latitude band")
            centers = vtk.vtkCellCenters(); centers.SetInputData(band); centers.Update()
            xy = np.array([centers.GetOutput().GetPoint(i)[:2]
                           for i in np.flatnonzero(active)])
            if reference_xy is None:
                reference_xy = xy
            elif xy.shape != reference_xy.shape or not np.allclose(xy, reference_xy, rtol=0, atol=1e-8):
                raise ValueError(f"{path}: layer columns do not align; remap to a common horizontal grid first")
            for key, name in (("dz", "dz"), ("Density", "Density"),
                              ("u", "Velocity_Zonal"), ("v", "Velocity_Meridional"),
                              ("w", "w"), ("Area", "Area")):
                fields[key].append(values(name)[active])
            if options.height_field:
                fields["z"].append(values(options.height_field)[active])
        heights = np.array(fields["z"]) if options.height_field else None
        result = energy_profile(*(np.array(fields[key]) for key in
                                  ("dz", "Density", "u", "v", "w", "Area")), options, heights)
        delta = result - mean
        mean += delta / count
        M2 += delta * (result - mean)
    std = np.sqrt(np.maximum(M2 / (count - 1), 0)) if count > 1 else np.zeros_like(mean)
    output.write("# EnergyDiagnostics: instantaneous area-weighted band anomalies, then time statistics.\n")
    output.write(f"# g={options.g}; rho_0={options.rho_0}; pressure_density={options.pressure_density}; surface_pressure={options.surface_pressure} Pa\n")
    output.write("# Pressure integrated upward through dz, shifted to specified surface gauge pressure.\n")
    output.write("# Complete water column required; cell-centred velocities assumed. No OMEGA fallback.\n")
    output.write("# Weights: polygon areas of supplied lon-lat mesh, as in existing diagnostics.\n")
    output.write("# ProfileConvergence is a derivative versus mean layer height, not a conservative moving-layer budget.\n")
    output.write("# Three-point nonuniform derivatives, one-sided endpoints; two layers use secants.\n")
    output.write("# Height: " + (options.height_field or "reconstructed from dz below flat z=0 surface") + "\n")
    output.write("# Pressure reference mode excludes density-anomaly contributions to hydrostatic pressure.\n")
    output.write("# Pressure fluxes use reduced p_d=p+rho_0*g*z; density anomalies retained by default.\n")
    output.write("# If surface elevation varies, provide --height-field with physical z to retain free-surface pressure effects.\n")
    output.write("# These five terms describe pressure transport, buoyancy conversion and vertical EKE advection, not a complete EKE budget.\n")
    energy_indices = list(range(1, len(ENERGY_TERMS)))
    height_index = 0
    output.write("# <q> is an instantaneous area-weighted band mean; q'=q-<q>.\n")
    output.write("# b=-g*(Density-rho_0)/rho_0; p_d=p+rho_0*g*z; z is positive upward.\n")
    output.write("# time_std is sample standard deviation across snapshots; zero for one snapshot.\n")
    output.write("# HeightCentre_time_mean_m: time mean of area-weighted centre height, in metres.\n")
    for i in energy_indices:
        name, units, definition = ENERGY_TERMS[i]
        output.write(f"# {name} [{units}]: {definition}\n")
    output.write("# layer HeightCentre_time_mean_m " + " ".join(
        ENERGY_TERMS[i][0] + suffix for i in energy_indices
        for suffix in ("_time_mean", "_time_std")) + "\n")
    for index in np.argsort(layers):
        numbers = np.column_stack((mean[index, energy_indices],
                                   std[index, energy_indices])).ravel()
        output.write(f"{layers[index]} {mean[index, height_index]:.10e} " +
                     " ".join(f"{x:.10e}" for x in numbers) + "\n")


def extract_band_with_area(data, lat1, lat2):
    """Use the same latitude selection and polygon areas as other statistics."""
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

    return band_with_area


def vertical_profile_derivative(flux, heights):
    """Differentiate in physical z, allowing either layer ordering."""
    flux = np.asarray(flux, dtype=float)
    heights = np.asarray(heights, dtype=float)
    if flux.ndim != 1 or heights.shape != flux.shape or flux.size < 2:
        raise ValueError("Vertical profile derivative needs at least two layer fluxes and heights")
    if not (np.all(np.isfinite(flux)) and np.all(np.isfinite(heights))):
        raise ValueError("Non-finite flux profile or layer heights")
    spacing = np.diff(heights)
    if not (np.all(spacing > 0.0) or np.all(spacing < 0.0)):
        raise ValueError("Mean layer heights must be strictly monotonic and distinct")
    return np.gradient(flux, heights, edge_order=2 if flux.size >= 3 else 1)


def compute_avg (data, field, lat1, lat2) :
    # Integrated zonal average statistics between latitudes lat1 and lat2

    band_with_area = extract_band_with_area(data, lat1, lat2)

    area_arr = band_with_area.GetCellData().GetArray("Area")

    if field_type == "Scalar":
        scalar = band_with_area.GetCellData().GetArray(field)
    
    zonal      = band_with_area.GetCellData().GetArray("Velocity_Zonal")
    meridional = band_with_area.GetCellData().GetArray("Velocity_Meridional")
    vertical   = get_vertical_velocity(band_with_area.GetCellData())
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
    vertical   = get_vertical_velocity(band_with_area.GetCellData())
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
if (len(sys.argv)<13) :
    print("Computes either averages or rms for given layers at given times. Returns mean and standard deviation over all times.")
    print("\nUsage: python statistics.py avg_type base_vtk_file k1 k2 t1 t2 dt lat1 lat2 scl_Omega scl_radius field_type [field | EnergyDiagnostics options]\n")
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
    print("field_type = Scalar, Rossby, KE, DeltaSM, DeltaI, VertFluxKE, EnergyDiagnostics")
    print("EnergyDiagnostics requires avg and a complete water column with at least two layers.")
    print("Layer numbering defaults to bottom-to-top; use --layer-order top-to-bottom to override.")
    print("Options: --g 9.80616 --rho-0 1030 --pressure-density in-situ --surface-pressure 0")
    print("Optional: --height-field Geopot_Height (physical centre height in metres).")
    print("field      = Field to analyze: \n \
               Options = \n \
                  Level \n \
                  Topography \n \
                  Penalization \n \
                  Ps/Eta \n \
                  Temperature/Density \n \
                  Velocity_Zonal \n \
                  Velocity_Meridional \n \
                  Vertical_Velocity/Velocity_Vertical \n \
                  Vorticity \n \
                  Geopot_Height \n \
                  P/Ps \n \
                  dz \n")
    print("Output is saved to run_[field_type]_[stat_type].txt")
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

valid_field_types = ("Scalar", "Rossby", "KE", "DeltaSM", "DeltaI", "VertFluxKE", "EnergyDiagnostics")
if field_type not in valid_field_types:
    sys.exit("ERROR: unsupported field_type. Choose " + ", ".join(valid_field_types))

if field_type == "EnergyDiagnostics":
    parser = argparse.ArgumentParser(description="Vertical energy fluxes and conversions with diagnosed hydrostatic pressure")
    parser.add_argument("--g", type=float, default=9.80616)
    parser.add_argument("--rho-0", type=float, default=1030.0)
    parser.add_argument("--pressure-density", choices=("reference", "in-situ"), default="in-situ")
    parser.add_argument("--surface-pressure", type=float, default=0.0, help="constant surface gauge pressure in Pa")
    parser.add_argument("--layer-order", choices=("top-to-bottom", "bottom-to-top"), default="bottom-to-top")
    parser.add_argument("--height-field", default=None, help="optional cell-centre physical z in m; otherwise reconstruct from dz below z=0")
    energy_options = parser.parse_args(sys.argv[13:])
    if stat_type != "avg" or k2 <= k1 or t2 < t1:
        sys.exit("ERROR: EnergyDiagnostics requires avg, at least two layers, and t2 >= t1")
    if not all(math.isfinite(x) for x in (energy_options.g, energy_options.rho_0, energy_options.surface_pressure)) or min(energy_options.g, energy_options.rho_0) <= 0:
        sys.exit("ERROR: g and rho_0 must be finite and positive; surface pressure must be finite")
    print("EnergyDiagnostics assumes the requested layers span the complete water column.", flush=True)

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

if field_type == "EnergyDiagnostics":
    try:
        write_energy_diagnostics(f, run, k1, k2, t1, t2, lat1, lat2, energy_options)
    except (ValueError, FileNotFoundError) as error:
        f.close()
        sys.exit(f"ERROR: {error}")

else:
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

f.close()

# Clean up temporary .vtk files (once, at the end)
for file in glob.glob(run+'_tri_[0-9][0-9][0-9]_[0-9][0-9][0-9][0-9].vtk'):
    os.remove(file)


