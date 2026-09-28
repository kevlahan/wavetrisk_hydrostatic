# Usage example (serial): python lonlat_to_3D.py SimpleJ5Z30 y 5 5 30 y fall 163 168
# Usage example (mpi):    srun -n 30 python -m mpi4py lonlat_to_3D.py SimpleJ5Z30 y 5 5 30 y fall 163 168 --step-days 5
#                         srun -n 30 python -m mpi4py lonlat_to_3D.py SimpleJ5Z30 y 5 5 30 n 163 168 --step-days 5
#    
#    Generates a 3D and zonal/meridional projections from a series of vtk layers.
#
#    Data has dimensions N x N/2 x NZ, where NZ is the number of vertical layers.
#    The vertical coordinate is P/Ps.
#
#    (1) vtkPoints are a concatenated array of vtkPoints
#    (2) vtkCell are constructed by joining the two 2D cells in adjacent layers to produce a vtkWedge prism.
#    (3) Each wedge retains its source layer's cell fields; interface pressure is
#        the average of point pressures from the two adjacent layers.
#
# Author: Weiguang Guan and Nicholas Kevlahan (McMaster University)
# Date : Last revision 2025-12-17 (Nicholas Kevlahan)

import csv, os, sys, glob, re, argparse, textwrap
from utilities import *
from contextlib import suppress
import numpy as np
import vtk
import subprocess
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtk.util import numpy_support
import scipy.ndimage
from datetime import datetime, timedelta
from pathlib import Path
from vtk.util.numpy_support import numpy_to_vtkIdTypeArray
from xyz2lonlat_layers import run_xyz_layers, run_xyz_layers_serial
from mpi4py import MPI

time_average_mode = "memory"
image_time_means = {}


# Missing physical samples are NaN; missing categorical Level is -1.
# Projected Level is the finest (maximum) valid contributing level.
def valid_mean(values, axis):
    good = np.isfinite(values)
    total = np.where(good, values, 0.0).sum(axis=axis, dtype=np.float64)
    count = good.sum(axis=axis)
    return np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)


def project_level(values, axis):
    return np.max(values, axis=axis).astype(np.int32)


def masked_vertical_smooth(values, valid, weights=None):
    values = np.asarray(values, dtype=np.float64)
    if valid.all():
        return scipy.ndimage.gaussian_filter1d(values, sigma=1, axis=0)
    if weights is None:
        weights = scipy.ndimage.gaussian_filter1d(valid.astype(float), sigma=1, axis=0)
    total = scipy.ndimage.gaussian_filter1d(np.where(valid, values, 0.0), sigma=1, axis=0)
    result = np.divide(total, weights, out=np.full_like(total, np.nan), where=weights > 0)
    result[~valid] = np.nan
    return result



def check_layer_correspondence(lower, upper):
    """Wedges require identical horizontal point order and triangle connectivity."""
    if (lower.GetNumberOfPoints() != upper.GetNumberOfPoints() or
            lower.GetNumberOfCells() != upper.GetNumberOfCells()):
        raise RuntimeError("Projected layers have different topology; regenerate all layers together")
    a = vtk_to_numpy(lower.GetPoints().GetData())[:, :2]
    b = vtk_to_numpy(upper.GetPoints().GetData())[:, :2]
    if (not np.allclose(a, b, rtol=0, atol=1e-8) or
            not np.array_equal(vtk_to_numpy(lower.GetPolys().GetData()),
                               vtk_to_numpy(upper.GetPolys().GetData()))):
        raise RuntimeError("Projected layer points/triangles do not correspond; cannot construct wedges")


################################################################################
class Cell3D():
    def __init__(self, vtp_series):
        global data_names
        self.vtp_series = vtp_series
        self.first = read_projected_layer(vtp_series[0])
        self.num_attrs = self.first.GetCellData().GetNumberOfArrays()
        self.attr_names = [self.first.GetCellData().GetArrayName(i) for i in range(self.num_attrs)]
        data_names = self.attr_names
        self.ugrid = vtk.vtkUnstructuredGrid()

    def construct(self):
        # Read each layer once. Allocate complete geometry and fields once.
        first = self.first
        npoint, ncell = first.GetNumberOfPoints(), first.GetNumberOfCells()
        nlayer = len(self.vtp_series)
        triangles = triangle_connectivity(first)
        xy = vtk_to_numpy(first.GetPoints().GetData())[:, :2]
        coordinates = np.empty(((nlayer + 1) * npoint, 3), dtype=np.float64)
        interfaces = coordinates.reshape(nlayer + 1, npoint, 3)
        interfaces[:, :, :2] = xy
        interfaces[0, :, 2] = 1.0
        interfaces[-1, :, 2] = 0.0
        connectivity = np.empty((nlayer * ncell, 6), dtype=np.int64)
        arrays = {}
        for name in self.attr_names:
            source = vtk_to_numpy(first.GetCellData().GetArray(name))
            arrays[name] = np.empty((nlayer * ncell,) + source.shape[1:], dtype=source.dtype)
        degree = np.bincount(triangles.ravel(), minlength=npoint)
        if np.any(degree == 0):
            raise RuntimeError("Projected mesh contains unused points")
        previous_pressure = None
        for k, path in enumerate(self.vtp_series):
            layer = first if k == 0 else read_projected_layer(path)
            if k:
                check_layer_correspondence(first, layer)
            pressure = point_cell_mean(layer, 'P/Ps', triangles, degree) * vert_scale
            if k:
                interfaces[k, :, 2] = (previous_pressure + pressure) * 0.5
            previous_pressure = pressure
            block = slice(k*ncell, (k+1)*ncell)
            connectivity[block, :3] = triangles + k*npoint
            connectivity[block, 3:] = triangles + (k+1)*npoint
            for name in self.attr_names:
                source_array = layer.GetCellData().GetArray(name)
                if source_array is None:
                    raise RuntimeError(f"Missing cell field {name!r} in {path}")
                arrays[name][block] = vtk_to_numpy(source_array)
        points = vtk.vtkPoints()
        points.SetData(numpy_to_vtk(coordinates, deep=False))
        self.ugrid.SetPoints(points)
        offsets = np.arange(0, connectivity.size + 1, 6, dtype=np.int64)
        cells = vtk.vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(offsets, deep=False),
                      numpy_to_vtkIdTypeArray(connectivity.ravel(), deep=False))
        self.ugrid.SetCells(vtk.VTK_WEDGE, cells)
        for name, values in arrays.items():
            array = numpy_to_vtk(values, deep=False)
            array.SetName(name)
            self.ugrid.GetCellData().AddArray(array)
        self.ugrid._numpy_buffers = [coordinates, connectivity, offsets, arrays]
        self.first = None
        delete_files(run+'_tri'+'*'+str(t).zfill(4)+'.vtp')
        return self.ugrid

    def construct_3Dimage(self):
        global Ntot
        global covarAvT, covarAvU, covarAvV, covarAvUV, covarAvVT
        global meanAvRho, meanAvT, meanAvU, meanAvV

        # Construct an unstructured grid
        ugrid = self.construct()

         # Write out 3D unstructured data
        if (len(idxs)==1):
            writer = vtk.vtkUnstructuredGridWriter()
            writer.SetFileTypeToBinary()
            writer.SetFileName(f"{run}_{t:04d}.vtk")
            writer.SetInputData(ugrid)
            writer.Write()

        # Resample ugrid to a regular grid
        img = vtk.vtkResampleToImage()
        img.SetInputDataObject(ugrid)        
        img.SetUseInputBounds(False)
        img.SetSamplingDimensions(lon_dim, lat_dim, vert_dim)
        img.SetSamplingBounds(lon_min, lon_max, lat_min, lat_max, vert_min, vert_max)
        img.Update()
        
        rgrid    = img.GetOutput()
        pnt_data = rgrid.GetPointData()

        img1 = vtk.vtkImageData()
        img2 = vtk.vtkImageData()
        img3 = vtk.vtkImageData()

        img1.SetDimensions(1,       lat_dim, vert_dim);
        img2.SetDimensions(lon_dim, 1,       vert_dim);
        img3.SetDimensions(lon_dim, lat_dim, vert_dim);

        img1.SetSpacing(rgrid.GetSpacing())
        img2.SetSpacing(rgrid.GetSpacing())
        img3.SetSpacing(rgrid.GetSpacing())

        img1.SetOrigin(rgrid.GetOrigin())
        img3.SetOrigin(rgrid.GetOrigin())

        shape = (vert_dim, lat_dim, lon_dim)
        mask_array = pnt_data.GetArray("vtkValidPointMask")
        if mask_array is None:
            raise RuntimeError("Resampler did not provide vtkValidPointMask")
        valid = vtk_to_numpy(mask_array).reshape(shape).astype(bool)
        fields = {}
        for name in self.attr_names:
            array = pnt_data.GetArray(name)
            if array is None or array.GetNumberOfComponents() != 1:
                raise RuntimeError(f"Missing or nonscalar input field: {name}")
            fields[name] = vtk_to_numpy(array).reshape(shape)
            valid &= np.isfinite(fields[name])
        if compressible == 'y':
            valid &= fields["Temperature"] > 0
        if "Level" in fields:
            levels = fields["Level"][valid]
            if np.any(levels < 0) or not np.allclose(levels, np.rint(levels), rtol=0, atol=1e-6):
                raise RuntimeError("Valid resampled Level values must be nonnegative integers")
        sample_counts = valid.sum(axis=2)
        for image, mask in ((img1, valid.any(axis=2)), (img2, valid.any(axis=1)), (img3, valid)):
            arr = numpy_to_vtk(mask.astype(np.uint8).ravel(), deep=True)
            arr.SetName("vtkValidPointMask")
            image.GetPointData().AddArray(arr)
        img2.SetOrigin(rgrid.GetOrigin())
        smooth_weights = None if valid.all() else scipy.ndimage.gaussian_filter1d(valid.astype(float), sigma=1, axis=0)
        vertical_profile = []
        for name in self.attr_names:
            if name == "Level":
                globe = np.where(valid, np.rint(fields[name]), -1).astype(np.int32)
                zonal = project_level(globe, 2)
                merid = project_level(globe, 1)
                profl = project_level(globe, (1, 2)).tolist()
            else:
                globe = masked_vertical_smooth(fields[name], valid, smooth_weights)
                zonal = valid_mean(globe, 2)
                merid = valid_mean(globe, 1)
                profl = valid_mean(globe, (1, 2)).tolist()
            for image, values in ((img1, zonal), (img2, merid), (img3, globe)):
                arr = numpy_to_vtk(values.ravel(), deep=True)
                arr.SetName(name)
                image.GetPointData().AddArray(arr)
            vertical_profile.append(profl)

        pnt_data = img3.GetPointData() # use smoothed data

        if compressible=='y':
            # Compute zonal projection of density and update average density
            meanRho   = compute_mean_density (pnt_data)
            meanAvRho = merge_mean(meanRho, sample_counts, meanAvRho, Ntot)

            # Compute covariances
            covarT, meanT, _ = compute_covar(pnt_data, "Temperature",         "Temperature")         # temperature 
            covarU, meanU, _ = compute_covar(pnt_data, "Velocity_Zonal",      "Velocity_Zonal")      # zonal velocity
            covarV, meanV, _ = compute_covar(pnt_data, "Velocity_Meridional", "Velocity_Meridional") # meridional velocity
            covarUV, _,    _ = compute_covar(pnt_data, "Velocity_Meridional", "Velocity_Zonal")      # momentum flux
            covarVT, _,    _ = compute_covar(pnt_data, "Velocity_Meridional", "Temperature")         # eddy heat flux
        
            # Update average covariances
            covarAvT,  newAvT, _ = merge_covariance(covarT,  meanT, meanT, sample_counts, covarAvT,  meanAvT, meanAvT, Ntot)
            covarAvU,  newAvU, _ = merge_covariance(covarU,  meanU, meanU, sample_counts, covarAvU,  meanAvU, meanAvU, Ntot)
            covarAvV,  newAvV, _ = merge_covariance(covarV,  meanV, meanV, sample_counts, covarAvV,  meanAvV, meanAvV, Ntot)
            covarAvUV, _,      _ = merge_covariance(covarUV, meanU, meanV, sample_counts, covarAvUV, meanAvU, meanAvV, Ntot)
            covarAvVT, _,      _ = merge_covariance(covarVT, meanV, meanT, sample_counts, covarAvVT, meanAvV, meanAvT, Ntot)

            # Update average means
            meanAvT = newAvT
            meanAvU = newAvU
            meanAvV = newAvV

            # Update total number of samples for averages
            Ntot    = Ntot + sample_counts
        
        file_out = f"{run}_{t:04d}"
        writer = vtk.vtkXMLImageDataWriter()
        images = (("", img3), ("_zonal", img1), ("_merid", img2))
        if time_average_mode == "memory":
            for suffix, image in images:
                if suffix not in image_time_means:
                    image_time_means[suffix] = ImageTimeMean(image)
                image_time_means[suffix].add(image)
        else:
            for suffix, image in images:
                writer.SetFileName(f"{file_out}{suffix}_{season}.vti")
                writer.SetInputData(image)
                if writer.Write() != 1:
                    raise RuntimeError("Failed to write snapshot image")

        if compressible=='y':
            if t == idxs[-1]: # write average statistics
                statistics = vtk.vtkImageData()
                statistics.SetDimensions(1, lat_dim, vert_dim);
                statistics.SetSpacing(rgrid.GetSpacing())
                statistics.SetOrigin(rgrid.GetOrigin())

                ke        = 0.5 * (covarAvU + covarAvV) # eddy kinetic energy
                mom_flux  = covarAvUV                   # eddy momentum flux
                heat_flux = covarAvVT                   # eddy heat flux                  

                add_scalar_data(np.where(Ntot > 0, covarAvT, np.nan), Nzonal, "TemperatureVariance", statistics)
                add_scalar_data(np.where(Ntot > 0, heat_flux, np.nan), Nzonal, "EddyHeatFlux",        statistics)
                add_scalar_data(np.where(Ntot > 0, mom_flux, np.nan), Nzonal, "EddyMomentumFlux",    statistics)
                add_scalar_data(np.where(Ntot > 0, ke, np.nan), Nzonal, "EddyKineticEnergy",   statistics)

                count_array = numpy_to_vtk(Ntot.ravel(), deep=True)
                count_array.SetName("ValidSampleCount")
                statistics.GetPointData().AddArray(count_array)
                mask = numpy_to_vtk((Ntot > 0).astype(np.uint8).ravel(), deep=True)
                mask.SetName("vtkValidPointMask")
                statistics.GetPointData().AddArray(mask)
                stats_file = f"{run}_statistics_zonal_{season}.vti"
                    
                writer.SetFileName(stats_file)
                writer.SetInputData(statistics)
                if writer.Write() != 1 or writer.GetErrorCode() != 0:
                    raise RuntimeError(f"Failed to write {stats_file}")

                # Save vertical profiles in csv file
                if len(idxs) == 1 or time_average_mode == "disk":
                    vertical_profile = np.array(vertical_profile).T
                    csv_file = f"{file_out}_{season}.csv"
                    with open(csv_file, 'w', newline='') as file:
                        writer = csv.writer(file)
                        writer.writerow(data_names)
                        for row in vertical_profile:
                            writer.writerow(row)

                
###############################################################################################################################
# Loads data from a vtk polydata file and sets vertical coordinate of interface based on attr_to_fill
def load_dataset(file1, file2, attr_to_fill=None):
    # Independent readers: reusing GetOutput() aliases the same VTK object.
    ugrid1 = read_projected_layer(file1)
    ugrid2 = read_projected_layer(file2)

    if (attr_to_fill is not None):
        set_vert_coord(ugrid1, ugrid2, attr_to_fill) # interface coordinate is average of coordinates of adjacent layers

    return ugrid2

#########################################################################################################################################
# Sets third coordinate of a vertex on specified interface as  average of data defined by attr_name of all cells that share the vertex.
# Currently, the vertical dimension is the scaled value of the average normalized pressure, P/Ps
def read_projected_layer(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    result = reader.GetOutput()
    if result.GetNumberOfCells() == 0:
        raise RuntimeError(f"No projected triangles in {path}")
    return result


def triangle_connectivity(layer):
    polys = layer.GetPolys()
    if not np.all(np.diff(vtk_to_numpy(polys.GetOffsetsArray())) == 3):
        raise RuntimeError("Projected layers must contain triangles")
    return vtk_to_numpy(polys.GetConnectivityArray()).reshape(-1, 3)


def point_cell_mean(layer, name, triangles=None, degree=None):
    if triangles is None:
        triangles = triangle_connectivity(layer)
    npoint = layer.GetNumberOfPoints()
    if degree is None:
        degree = np.bincount(triangles.ravel(), minlength=npoint)
    values = vtk_to_numpy(layer.GetCellData().GetArray(name)).astype(np.float64, copy=False)
    total = np.zeros(npoint, dtype=np.float64)
    for vertex in range(3):
        total += np.bincount(triangles[:, vertex], weights=values, minlength=npoint)
    if np.any(degree == 0):
        raise RuntimeError("Projected mesh contains unused points")
    return total / degree


def set_vert_coord(ugrid1, ugrid2, attr_name):
    check_layer_correspondence(ugrid1, ugrid2)
    coords = vtk_to_numpy(ugrid1.GetPoints().GetData()).astype(np.float64, copy=True)
    coords[:, 2] = (point_cell_mean(ugrid1, attr_name) + point_cell_mean(ugrid2, attr_name)) * (vert_scale / 2)
    ugrid2.GetPoints().SetData(numpy_to_vtk(coords, deep=True))


# Loads data from a vtk polydata file and sets vertical coordinate of top and bottom interfaces based on attr_to_fill
def load_dataset_bdry(file, interface, attr_to_fill=None):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file)
    reader.Update()
    ugrid = reader.GetOutput()

    if (attr_to_fill is not None):
        set_vert_coord_bdry(ugrid, interface) # set vertical coordinate to be scaled value of data given by attr_to_fill (P/Ps)

    return ugrid


#########################################################################################################################################
# Sets third coordinate of a vertex for top and bottom interfaces
def set_vert_coord_bdry(ugrid, interface):
    # Get point coordinate array
    coords = vtk_to_numpy(ugrid.GetPoints().GetData())

    if interface == 0:    # ground
        coords[:,2] = 1 
    elif interface == nz: # top of atmosphere
        coords[:,2] = 0 

    # Update vertical coordinate for this layer
    ugrid.GetPoints().SetData(numpy_to_vtk(coords))
    

#########################################################################################################################################
# Functions to convert between images and numpy arrays
def read_vti_images(file_list):
    """
    Reads VTI image files from a list of file names.
    
    Args:
        file_list (list): List of VTI file paths to read.
    
    Returns:
        list: A list of VTK image data objects.
    """
    vti_images = []

    for file_name in file_list:
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(file_name)
        reader.Update()
        vti_images.append(reader.GetOutput())

    return vti_images


def time_mean():
    # Computes time means of meridional and zonal averages and vertical profile
    zonal = avg_images("_zonal") # zonal projection
    merid = avg_images("_merid") # meridional projection
    avg_images("")               # complete 3D data

    dims       = zonal.GetDimensions() 
    point_data = zonal.GetPointData()
    num_attrs  = point_data.GetNumberOfArrays()

    data_names       = []
    vertical_profile = []
    for i in range(num_attrs):
        data_names.append(point_data.GetArrayName(i))
        array = vtk_to_numpy(point_data.GetArray(i))
        array = array.reshape((vert_dim, lat_dim))
        profile = (project_level(array, 1) if point_data.GetArrayName(i) == "Level" else valid_mean(array, 1)).tolist()
        vertical_profile.append(profile)

    vertical_profile = np.array(vertical_profile).T
    
    with open(run+"_"+season+".csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data_names)
        for row in vertical_profile:
            writer.writerow(row)
   
    # Remove numbered VTI snapshots, including leftovers from previous runs,
    # only after all final averages and the profile have been written.
    remove_intermediate_vti()

    # Existing cleanup for other temporary outputs in disk mode.
    if len(idxs) != 1 and time_average_mode == "disk":
        for i in idxs:
             tag = f"_{i:04d}"
             delete_files(f"{run}{tag}_{season}.vtk")
             delete_files(f"{run}{tag}_{season}.csv")
             delete_files(f"{run}{tag}_zonal_{season}.vti")
             delete_files(f"{run}{tag}_merid_{season}.vti")
             delete_files(f"{run}{tag}_{season}.vti")


def remove_intermediate_vti():
    prefix = Path(run)
    pattern = re.compile(
        re.escape(prefix.name) + r"_\d+(?:_zonal|_merid|_statistics_zonal)?_"
        + re.escape(season) + r"\.vti$"
    )
    for path in prefix.parent.iterdir():
        if path.is_file() and pattern.fullmatch(path.name):
            path.unlink()

class ImageTimeMean:
    def __init__(self, image):
        self.geometry = vtk.vtkImageData()
        self.geometry.CopyStructure(image)
        pd = image.GetPointData()
        self.names = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())
                      if pd.GetArrayName(i) != "vtkValidPointMask"]
        n = image.GetNumberOfPoints()
        self.sums = {name: np.zeros(n, dtype=np.float64) for name in self.names if name != "Level"}
        self.counts = {name: np.zeros(n, dtype=np.uint32) for name in self.sums}
        self.finest = np.full(n, -1, dtype=np.int32) if "Level" in self.names else None
        self.valid = np.zeros(n, dtype=bool)

    def add(self, image):
        if (image.GetExtent() != self.geometry.GetExtent() or
                image.GetOrigin() != self.geometry.GetOrigin() or
                image.GetSpacing() != self.geometry.GetSpacing()):
            raise RuntimeError("Cannot average images on different grids")
        pd = image.GetPointData()
        mask = vtk_to_numpy(pd.GetArray("vtkValidPointMask")).astype(bool)
        self.valid |= mask
        for name in self.names:
            values = vtk_to_numpy(pd.GetArray(name))
            good = mask & np.isfinite(values)
            if name == "Level":
                np.maximum(self.finest, np.where(good, values, -1).astype(np.int32), out=self.finest)
            else:
                self.sums[name] += np.where(good, values, 0)
                self.counts[name] += good

    def finish(self):
        out = vtk.vtkImageData()
        out.CopyStructure(self.geometry)
        for name in self.names:
            if name == "Level":
                result = self.finest
            else:
                result = self.sums.pop(name)
                count = self.counts.pop(name)
                np.divide(result, count, out=result, where=count > 0)
                result[count == 0] = np.nan
            arr = numpy_to_vtk(result, deep=True)
            arr.SetName(name)
            out.GetPointData().AddArray(arr)
        mask = numpy_to_vtk(self.valid.astype(np.uint8), deep=True)
        mask.SetName("vtkValidPointMask")
        out.GetPointData().AddArray(mask)
        return out


def avg_images(file_type):
    if time_average_mode == "memory":
        output_image = image_time_means.pop(file_type).finish()
    else:
        accumulator = None
        reader = vtk.vtkXMLImageDataReader()
        for t in idxs:
            reader.SetFileName(run+f"_{t:04d}"+file_type+"_"+season+".vti")
            reader.Update()
            image = reader.GetOutput()
            if accumulator is None:
                accumulator = ImageTimeMean(image)
            accumulator.add(image)
        output_image = accumulator.finish()

    # ---- Write the averaged image ----
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(run+file_type+'_'+season+'.vti')
    writer.SetInputData(output_image)
    if writer.Write() != 1 or writer.GetErrorCode() != 0:
        raise RuntimeError(f"Failed to write final average for {file_type!r}")

    return output_image
    
def mean (data):
    # Computes mean of data

    mean = n = 0
    for x, in zip(data):
        n    += 1
        mean += (x - mean) / n
        
    return mean

def covariance(data1, data2):
    # Population covariance of available longitude samples (climatological moment).
    good = np.isfinite(data1) & np.isfinite(data2)
    x, y = np.asarray(data1)[good], np.asarray(data2)[good]
    if x.size == 0:
        return 0.0, 0.0, 0.0
    mx, my = x.mean(), y.mean()
    return np.mean((x - mx) * (y - my)), mx, my


def merge_mean(mean1, n1, mean2, n2):
    n = n1 + n2
    total = np.where(n1 > 0, n1 * mean1, 0) + np.where(n2 > 0, n2 * mean2, 0)
    return np.divide(total, n, out=np.zeros_like(total, dtype=float), where=n > 0)


def merge_covariance(cov1, mean1_x, mean1_y, n1,
                     cov2, mean2_x, mean2_y, n2):
    # Merge population moments using the valid count at each latitude/level.
    n = n1 + n2
    dx, dy = mean1_x - mean2_x, mean1_y - mean2_y
    cross = np.divide(n1 * n2, n, out=np.zeros_like(dx, dtype=float), where=n > 0)
    total = n1 * cov1 + n2 * cov2 + cross * dx * dy
    cov = np.divide(total, n, out=np.zeros_like(total, dtype=float), where=n > 0)
    return cov, merge_mean(mean1_x, n1, mean2_x, n2), merge_mean(mean1_y, n1, mean2_y, n2)


def compute_covar(pnt_data, var1, var2):
    shape = (vert_dim, lat_dim, lon_dim)
    x = vtk_to_numpy(pnt_data.GetArray(var1)).reshape(shape)
    y = vtk_to_numpy(pnt_data.GetArray(var2)).reshape(shape)
    good = np.isfinite(x) & np.isfinite(y)
    count = good.sum(axis=2)
    mx = np.divide(np.where(good, x, 0).sum(axis=2), count,
                   out=np.zeros(shape[:2]), where=count > 0)
    my = mx if var1 == var2 else np.divide(np.where(good, y, 0).sum(axis=2), count,
                                         out=np.zeros(shape[:2]), where=count > 0)
    dx = np.where(good, x - mx[:, :, None], 0)
    dy = dx if var1 == var2 else np.where(good, y - my[:, :, None], 0)
    cov = np.divide(np.einsum('ijk,ijk->ij', dx, dy), count,
                    out=np.zeros(shape[:2]), where=count > 0)
    return cov, mx, my


def compute_mean_density(pnt_data):
    shape = (vert_dim, lat_dim, lon_dim)
    T = vtk_to_numpy(pnt_data.GetArray("Temperature")).reshape(shape)
    Ps = vtk_to_numpy(pnt_data.GetArray("Ps")).reshape(shape)
    P_Ps = vtk_to_numpy(pnt_data.GetArray("P/Ps")).reshape(shape)
    return valid_mean(P_Ps * Ps / (Rd * T), 2)


def add_scalar_data(data, N, name, img):
    # Adds scalar data of total size N to vtk img with given name
    
    attr = numpy_to_vtk(data.reshape(N))
    attr.SetName(name)
    
    img.GetPointData().AddArray(attr)

class _Fmt(argparse.ArgumentDefaultsHelpFormatter,
           argparse.RawDescriptionHelpFormatter):
    pass

DOC = textwrap.dedent("""\
  Example 1:
    python lonlat_to_3D.py SimpleJ5J7Z30 y 5 7 30 n 1 365
      processes compressible data from run SimpleJ5J7Z30 with levels 5 to 7 and 30 layers
      for *.vtk.tgz data files with indices 1 to 365

  Example 2 (seasonal statistics):
    python lonlat_to_3D.py SimpleJ5J7Z30 y 5 7 30 y spring 1 365 --step-days 5
      processes compressible data from run SimpleJ5J7Z30 with levels 5 to 7 and 30 layers for spring
      for *.vtk.tgz data with indices 1 to 365 saved every 5 days

  Example 3 (seasonal statistics starting from June 21):
    python lonlat_to_3D.py SimpleJ5J7Z30 y 5 7 30 y spring 1 365 --step-days 5 --start-date 0001-06-21
      processes compressible data from run SimpleJ5J7Z30 with levels 5 to 7 and 30 layers for spring
      for *.vtk.tgz data with indices 1 to 365 saved every 5 days with first data from June 21.

  Final VTI outputs (also for a single selected time):
    <run>_<season>.vti                   3D time mean
    <run>_zonal_<season>.vti             zonal time mean
    <run>_merid_<season>.vti             meridional time mean
    <run>_statistics_zonal_<season>.vti  statistics (compressible data)

  Numbered intermediate VTI files for this run and season are removed after
  successful final output, including leftovers from previous invocations.
  Other runs and seasons are preserved. With seasonal mode disabled,
  <season> is 'all'. Single-time runs also retain the unstructured .vtk output.

  3D data has dimensions N x N/2 x nz. The vertical coordinate is P/Ps.
""")

DEFAULT_START_DATE_STR = "0001-03-22"

def parse_start_date(s: str) -> datetime:
    """
    Accepts:
      - YYYY-MM-DD (ISO), e.g., 2025-03-22
      - M-D (year defaults to 0001), e.g., 3-22 or 03-22
      - M/D likewise, e.g., 3/22
    Returns a datetime at 00:00:00.
    """
    if not s:
        return datetime(1, 3, 22)
    # try ISO first
    try:
        return datetime.fromisoformat(s)
    except ValueError:
        pass
    # allow M-D or M/D with default year=0001
    s2 = s.replace("/", "-")
    try:
        m, d = map(int, s2.split("-"))
        return datetime(1, m, d)
    except Exception:
        raise argparse.ArgumentTypeError("Use YYYY-MM-DD or M-D (year defaults to 0001).")

def parse_args():
    global time_average_mode
    p = argparse.ArgumentParser(
        prog="lonlat_to_3D.py",
        description=
        """
Generates a 3D data files, zonal/meridional projections and vertical profiles
from a series of layers in directory folder.
        
        """,
        epilog=DOC,
        formatter_class=_Fmt,
    )
    
    # required positionals
    p.add_argument("run")
    p.add_argument("compressible", choices=["y","n","Y","N"])
    p.add_argument("Jmin", type=int)
    p.add_argument("Jmax", type=int)
    p.add_argument("nz",   type=int)

    # seasons flag (y/n), then conditionally the season name
    p.add_argument("seasons", choices=["y","n","Y","N"],
                   help="y => seasonal mode (needs season + step-days); n => all")
    p.add_argument("season", nargs="?", default="all",
                   help="season name (required if seasons=='y')")
    # common time bounds
    p.add_argument("t1", type=int)
    p.add_argument("t2", type=int)

    # only used in seasonal mode; optional flag is more robust than a positional
    p.add_argument("--step-days", type=int, default=5,
                   help="day stride (seasonal mode); default 5")

    p.add_argument(
        "--start-date",
        type=parse_start_date,
        default=parse_start_date(DEFAULT_START_DATE_STR),
        help="reference start date. Formats: YYYY-MM-DD or M-D (year defaults to 0001).",
    )
    
    if len(sys.argv) == 1:  # no args → show full help + examples
        print()
        p.print_help(); sys.exit(0)

    p.add_argument("--time-average-mode", choices=("memory", "disk"), default="memory",
                   help="memory avoids temporary VTI files; disk keeps the averaging pass separate to lower peak memory")
    a = p.parse_args()
    time_average_mode = a.time_average_mode
    
    seasons = a.seasons.lower() == "y"
    if seasons and (a.season == "all"):
        p.error("seasonal mode requires a season name (e.g., 'spring').")

    if rank==0:
        for k, v in sorted(vars(a).items()):
            print(f"{k:>12}: {v}", flush=True)

    # unpack to variable names
    run          = a.run
    compressible = a.compressible
    Jmin         = a.Jmin
    Jmax         = a.Jmax
    nz           = a.nz
    t1           = a.t1
    t2           = a.t2
    season       = a.season if seasons else "all"
    step_days    = a.step_days if seasons else 5
    start_date   = a.start_date 

    if seasons:
        total_days = t2 * step_days
        n_years    = total_days // 365
    else:
        total_days = None
        n_years    = None

    return run, compressible, Jmin, Jmax, nz, seasons, season, t1, t2, step_days, start_date, total_days, n_years

    
#########################################################################################################################################
#    Main program
#########################################################################################################################################
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Input parameters
(run, compressible, Jmin, Jmax, nz, seasons, season, t1, t2, step_days, start_date, total_days, n_years) = parse_args()

# Grid dimensions (same number of rectangular cells as lozenge cells on the sphere) 
lat_dim  = int(np.sqrt((10*4**Jmax + 2)/2))
lon_dim  = 2*lat_dim
vert_dim = nz

dtheta_min = 180/lat_dim
dtheta_max = dtheta_min * 2**(Jmax - Jmin)

lon_min  = -180 
lon_max  =  180 
lat_min  =  -90 
lat_max  =   90 
vert_min =    0 
vert_max =    1

Ngrid    = lat_dim * lon_dim * vert_dim                                 
Nmerid   = lon_dim * vert_dim
Nzonal   = lat_dim * vert_dim

# Physical constants
Rd          = 287  # gas constant
H           = 4000 # mean depth
ref_density = 1030 # incompressible case
vert_scale  = 1    # rescaling for P/Ps vertical coordinate

vert_min = vert_scale * vert_min
vert_max = vert_scale * vert_max

# Initialize statistics variables
Ntot      = np.zeros((vert_dim, lat_dim), dtype=np.int64)
covarAvT  = np.zeros((vert_dim,lat_dim))
covarAvU  = np.zeros((vert_dim,lat_dim))
covarAvV  = np.zeros((vert_dim,lat_dim))
covarAvUV = np.zeros((vert_dim,lat_dim))
covarAvVT = np.zeros((vert_dim,lat_dim))
meanAvRho = np.zeros((vert_dim,lat_dim))
meanAvT   = np.zeros((vert_dim,lat_dim))
meanAvU   = np.zeros((vert_dim,lat_dim))
meanAvV   = np.zeros((vert_dim,lat_dim))

# Remove .DS_store to avoid load error
with suppress(OSError):
    os.remove(sys.argv[1]+'/.DS_Store')

if rank == 0:
    print(f"\nInterpolating to uniform {lon_dim} x {lat_dim} x {vert_dim} grid\n", flush=True)
    print(f"Season = {season}\n", flush=True)

if seasons:
    def ceil_div(a, b):
        return -(-a // b)

    # Start dates only; end is the next season's start 
    SEASON_START = {
        'spring': (3, 1),   # Mar 1
        'summer': (6, 1),   # Jun 1
        'fall':   (9, 1),   # Sep 1
        'winter': (12, 1),  # Dec 1
    }
    
    SEASON_ORDER = ['spring', 'summer', 'fall', 'winter']
    
    SEASON_MONTHS = {
        'spring': {3, 4, 5},
        'summer': {6, 7, 8},
        'fall':   {9, 10, 11},
        'winter': {12, 1, 2},
    }

    if season not in SEASON_MONTHS:
        raise ValueError("season must be one of: spring, summer, fall, winter")

    allowed = SEASON_MONTHS[season]
    idxs = []
    for t in range(t1, t2 + 1):  
        dt = start_date + timedelta(days=t * step_days)
        if dt.month in allowed:
            idxs.append(t)
else:
    idxs = range (t1, t2 + 1)
    
idxs = sorted(set(idxs)) # remove duplicates and sort

if rank == 0:
     print(f"{len(idxs)} file indices to process: {list(idxs)}\n", flush=True)

for t in idxs:
    if rank == 0:
        print(f"    Processing file with index {t}", flush=True)
        
    failures = run_xyz_layers(run, Jmin, Jmax, nz, t, t, True)
    failed = comm.bcast(bool(failures) if rank == 0 else None, root=0)
    if failed:
        raise RuntimeError(f"Layer projection failed: {failures}")

    comm.Barrier()          # ensure all .vtp files are finished

    if rank != 0:
        continue            # only rank 0 continues past here
   
    vtp_series = []
    for z in range (1, nz+1):
        vtp_series.append(run+'_tri_lonlat_'+str(z).zfill(3)+'_'+str(t).zfill(4)+".vtp")
            
    cell3d = Cell3D(vtp_series)
    cell3d.construct_3Dimage()

# Compute mean over all times
if rank == 0 and idxs:
    time_mean()

#########################################################################################################################################
