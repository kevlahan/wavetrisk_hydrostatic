# Usage: python lonlat_to_3D.py run nz t1 t2 J 
#    
#    Generates a 3D and zonal/meridional projections from a series of vtk layers.
#
#    Data has dimensions N x N/2 x NZ, where NZ is the number of vertical layers.
#    The vertical coordinate is P/Ps.
#
#    (1) vtkPoints are a concatenated array of vtkPoints
#    (2) vtkCell are constructed by joining the two 2D cells in adjacent layers to produce a vtkWedge prism.
#    (3) vtkCellData are computed as the average of the values from the two 2D cells in adjacent layers.  
#
# Author: Weiguang Guan and Nicholas Kevlahan (McMaster University)
# Date  : Last revision 2025-01-30 (Nicholas Kevlahan)

import os
import sys
from utilities import *
from contextlib import suppress
import numpy as np
import vtk
import csv
import subprocess
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import scipy.ndimage

################################################################################
class Cell3D() :
    def __init__(self, vtp_series) :
        global data_names
        self.vtp_series = vtp_series
        self.ugrid = vtk.vtkUnstructuredGrid()  # final result

        # Load DS1
        ds1 = load_dataset_bdry(vtp_series[0], 0, 'P/Ps')

        points1 = ds1.GetPoints()
        self.ugrid.SetPoints(points1)

        self.num_attrs = ds1.GetCellData().GetNumberOfArrays()

        self.attr_list  = [None] * self.num_attrs     # to store attribute arrays
        self.attr_names = [None] * self.num_attrs

        for i in range(self.num_attrs) :
            attr1 = ds1.GetCellData().GetArray(i)
            self.attr_names[i] = attr1.GetName()

        data_names = self.attr_names
        
        
    # Progressively constructs 3D wedge cells, from surface interface upwards   
    def construct(self) :
        # Bottom interface is ground: P/Ps=1 by definition
        file1 = self.vtp_series[0]
        ds1   = load_dataset_bdry(file1, 0, 'P/Ps')
        
        for i, file2 in enumerate(self.vtp_series[1:]) : # loop through layers to (file2 is used to determine upper interface)
            ds2 = load_dataset(file1, file2, 'P/Ps') # upper interface

            # Construct cells between two layers and add them to ugrid
            self.add_cells(ds1, ds2, (i+1))

            ds1   = ds2 # set upper interface as next lower interfaces
            file1 = file2

        # Top interface is P/Ps=0 by definition
        ds2 = load_dataset_bdry(file2, nz, 'P/Ps')
        self.add_cells(ds1, ds2, (i+1))

        # Set attributes for the ugrid
        for i in range(self.num_attrs) :
            attr_array = numpy_to_vtk(self.attr_list[i])
            attr_array.SetName(self.attr_names[i])

            self.ugrid.GetCellData().AddArray(attr_array)
            
        vtp_files = run+'_tri'+"*"+str(t).zfill(4)+".vtp"
        delete_files(vtp_files)

        return self.ugrid

    # Function to construct cells between two layers (ds1 and ds2) and
    # then add them to the unstructured grid (ugrid)
    # Here, ugrid has already had vtkPoints from ds1
    def add_cells(self, ds1, ds2, cell_layer_id=0) :
        # Points and coordinates at all interfaces below current layer
        points         = self.ugrid.GetPoints()
        coords         = vtk_to_numpy(points.GetData())

        # Upper interface
        points2        = ds2.GetPoints()
        num_points2    = ds2.GetNumberOfPoints()
        coords2        = vtk_to_numpy(points2.GetData())
        num_cells2     = ds2.GetNumberOfCells()
        cellformation2 = vtk_to_numpy(ds2.GetPolys().GetData())
            
        # Add points and coordinates for upper interface
        coords = np.concatenate((coords, coords2))
        points.SetData(numpy_to_vtk(coords))
        self.ugrid.SetPoints(points)
        num_points  = points.GetNumberOfPoints()
        
        startID = 0
        # loop through cells in ds2
        for i in range(num_cells2) :
            size     = cellformation2[startID]
            pnt_ids  = cellformation2[(startID+1):(startID+1+size)]

            pnt_ids1 = pnt_ids + (num_points - 2*num_points2)
            pnt_ids2 = pnt_ids + (num_points -   num_points2)
            
            # 6 points in wedge cell
            pnt_ids  = np.concatenate((pnt_ids1, pnt_ids2))

            # Triangular cells
            self.ugrid.InsertNextCell(vtk.VTK_WEDGE , len(pnt_ids), pnt_ids)
    
            startID = startID + 1 + size

        # Set the cell attributes from layer between lower and upper interfaces
        for i in range(self.num_attrs) :
            attr = ds1.GetCellData().GetArray(i)

            if (self.attr_list[i] is None) :
                self.attr_list[i] = attr
            else :
                self.attr_list[i] = np.concatenate((self.attr_list[i], attr))
    
    def construct_3Dimage(self) :
        global Ntot
        global covarAvT, covarAvU, covarAvV, covarAvUV, covarAvVT
        global meanAvRho, meanAvT, meanAvU, meanAvV

        # Construct an unstructured grid
        ugrid = self.construct()

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

        img1.SetDimensions(1, lat_dim, vert_dim);
        img2.SetDimensions(lon_dim, 1, vert_dim);
        img3.SetDimensions(lon_dim, lat_dim, vert_dim);

        img1.SetSpacing(rgrid.GetSpacing())
        img2.SetSpacing(rgrid.GetSpacing())
        img3.SetSpacing(rgrid.GetSpacing())

        img1.SetOrigin(rgrid.GetOrigin())
        img3.SetOrigin(rgrid.GetOrigin())

        vertical_profile = []
        for i in range(self.num_attrs) :
            name = pnt_data.GetArray(i).GetName()

            # Extract 3D global data
            globe = vtk_to_numpy(pnt_data.GetArray(i))
            globe = globe.reshape(vert_dim, lat_dim, lon_dim)
            globe = scipy.ndimage.gaussian_filter(globe, sigma=(1.0, 0, 0)) 

            zonal = np.mean(globe, axis=2)              # zonal average
            merid = np.mean(globe, axis=1)              # meridional average
            profl = np.mean(globe, axis=(1,2)).tolist() # vertical profile

            attr1 = numpy_to_vtk(zonal.reshape(Nzonal)); attr1.SetName(name)
            attr2 = numpy_to_vtk(merid.reshape(Nmerid)); attr2.SetName(name)
            attr3 = numpy_to_vtk(globe.reshape(Ngrid));  attr3.SetName(name)

            img1.GetPointData().AddArray(attr1)
            img2.GetPointData().AddArray(attr2)
            img3.GetPointData().AddArray(attr3)
            vertical_profile.append(profl) 

        # Write out 3D unstructured data
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
        writer.SetFileName(sys.argv[1]+"_"+str(t).zfill(4)+".vtk")
        writer.SetInputData(ugrid)
        writer.Write()
        
        pnt_data = img3.GetPointData() # use smoothed data

        if compressible=='y' :
            # Compute zonal projection of density and update average density
            meanRho   = compute_mean_density (pnt_data)
            meanAvRho = merge_mean(meanRho, lon_dim, meanAvRho, Ntot)

            # Compute covariances
            covarT, meanT, _ = compute_covar(pnt_data, "Temperature",        "Temperature")        # temperature 
            covarU, meanU, _ = compute_covar(pnt_data, "VelocityZonal",      "VelocityZonal")      # zonal velocity
            covarV, meanV, _ = compute_covar(pnt_data, "VelocityMeridional", "VelocityMeridional") # meridional velcoity
            covarUV, _,    _ = compute_covar(pnt_data, "VelocityZonal",      "VelocityMeridional") # momentum flux
            covarVT, _,    _ = compute_covar(pnt_data, "VelocityMeridional", "Temperature")        # eddy heat flux
        
            # Update average covariances
            covarAvT,  newAvT, _ = merge_covariance(covarT,  meanT, meanT, lon_dim, covarAvT,  meanAvT, meanAvT, Ntot)
            covarAvU,  newAvU, _ = merge_covariance(covarU,  meanU, meanU, lon_dim, covarAvU,  meanAvU, meanAvU, Ntot)
            covarAvV,  newAvV, _ = merge_covariance(covarV,  meanV, meanV, lon_dim, covarAvV,  meanAvV, meanAvV, Ntot)
            covarAvUV, _,      _ = merge_covariance(covarUV, meanU, meanV, lon_dim, covarAvUV, meanAvU, meanAvV, Ntot)
            covarAvVT, _,      _ = merge_covariance(covarVT, meanV, meanT, lon_dim, covarAvVT, meanAvV, meanAvT, Ntot)

            # Update average means
            meanAvT = newAvT
            meanAvU = newAvU
            meanAvV = newAvV

            # Update total number of samples for averages
            Ntot    = Ntot + lon_dim
        
        # Write out data
        file_out = run+"_"+str(t).zfill(4)

        # Write image (uniform grid) data
        writer = vtk.vtkXMLImageDataWriter()

        # Write 3D Cartesian grid data
        writer.SetFileName(file_out+".vti")
        writer.SetInputData(img3)
        writer.Write()

        # Write zonal projection
        writer.SetFileName(file_out+"_zonal.vti")
        writer.SetInputData(img1)
        writer.Write()

        # Write meridional projection
        writer.SetFileName(file_out+"_merid.vti")
        writer.SetInputData(img2)
        writer.Write()

        if compressible=='y' :
            if t == t2: # write average statistics
                statistics = vtk.vtkImageData()
                statistics.SetDimensions(1, lat_dim, vert_dim);
                statistics.SetSpacing(rgrid.GetSpacing())
                statistics.SetOrigin(rgrid.GetOrigin())

                ke        = 0.5 * meanAvRho * (covarAvU + covarAvV) # eddy kinetic energy
                mom_flux  = 0.5 * meanAvRho * covarAvUV             # eddy momentum flux
                heat_flux =       meanAvRho * covarAvVT             # eddy heat flux                  

                add_scalar_data(covarAvT,  Nzonal, "TemperatureVariance", statistics)
                add_scalar_data(heat_flux, Nzonal, "EddyHeatFlux",        statistics)
                add_scalar_data(mom_flux,  Nzonal, "EddyMomentumFlux",    statistics)
                add_scalar_data(ke,        Nzonal, "EddyKineticEnergy",   statistics)
                
                writer.SetFileName(run+"_statistics_mean.vti")
                writer.SetInputData(statistics)
                writer.Write()

                # Save vertical profiles in csv file
                vertical_profile = np.array(vertical_profile).T
                with open(file_out+"_profile.csv", 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(data_names)
                    for row in vertical_profile :
                        writer.writerow(row)

                
###############################################################################################################################
# Loads data from a vtk polydata file and sets vertical coordinate of interface based on attr_to_fill
def load_dataset(file1, file2, attr_to_fill=None) :
    reader = vtk.vtkXMLPolyDataReader()
        
    # Layer below current interface
    reader.SetFileName(file1)
    reader.Update()
    ugrid1 = reader.GetOutput() 

    # Layer above current interface
    reader.SetFileName(file2)
    reader.Update()
    ugrid2 = reader.GetOutput()

    if (attr_to_fill is not None) :
        set_vert_coord(ugrid1, ugrid2, attr_to_fill) # interface coordinate is average of coordinates of adjacent layers

    return ugrid2

#########################################################################################################################################
# Sets third coordinate of a vertex on specified interface as  average of data defined by attr_name of all cells that share the vertex.
# Currently, the vertical dimension is the scaled value of the average normalized pressure, P/Ps
def set_vert_coord(ugrid1, ugrid2, attr_name) :
    # Load data from layers adjacent to current interface
    coords1        = vtk_to_numpy(ugrid1.GetPoints().GetData())
    num_points1    = coords1.shape[0]
    attrs1         = vtk_to_numpy(ugrid1.GetCellData().GetArray(attr_name)) # data to use as the vertical coordinate (currently P/Ps)
    cellformation1 = vtk_to_numpy(ugrid1.GetPolys().GetData())
    num_cells1     = ugrid1.GetNumberOfCells()
        
    coords2        = vtk_to_numpy(ugrid2.GetPoints().GetData())
    num_points2    = coords2.shape[0]
    attrs2         = vtk_to_numpy(ugrid2.GetCellData().GetArray(attr_name)) 
    cellformation2 = vtk_to_numpy(ugrid2.GetPolys().GetData())
    num_cells2     = ugrid2.GetNumberOfCells()   

    # Loop through in layer below interface
    weights = np.zeros(num_points1) 
    startID = 0
    for i in range(num_cells1) :
        size    = cellformation1[startID]
        pnt_ids = cellformation1[(startID+1):(startID+1+size)]
        for pnt in pnt_ids :
            coords1[pnt,2]  = coords1[pnt,2] * weights[pnt] + attrs1[i] * vert_scale # compute average value of data used for vertical coordinate
            weights[pnt]   += 1
            coords1[pnt,2] /= weights[pnt]                                           # increment weighted average
        startID = startID + 1 + size

    # Loop through cells in layer above interfaces
    startID = 0
    weights = np.zeros(num_points2) 
    for i in range(num_cells2) :
        size    = cellformation2[startID]
        pnt_ids = cellformation2[(startID+1):(startID+1+size)]
        for pnt in pnt_ids :
            coords2[pnt,2] = coords2[pnt,2] * weights[pnt] + attrs2[i] * vert_scale # compute average value of data used for vertical coordinate
            weights[pnt]   += 1
            coords2[pnt,2] /= weights[pnt]                                          # increment weighted average
        startID = startID + 1 + size

    # Interface vertical coordinate is average of vertical coordinates of adjacent layers
    coords1[:,2] = (coords1[:,2] + coords2[:,2])/2.0
    
    # Update vertical coordinate for this layer
    ugrid2.GetPoints().SetData(numpy_to_vtk(coords1))


# Loads data from a vtk polydata file and sets vertical coordinate of top and bottom interfaces based on attr_to_fill
def load_dataset_bdry(file, interface, attr_to_fill=None) :
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file)
    reader.Update()
    ugrid = reader.GetOutput()

    if (attr_to_fill is not None) :
        set_vert_coord_bdry(ugrid, interface) # set vertical coordinate to be scaled value of data given by attr_to_fill (P/Ps)

    return ugrid


#########################################################################################################################################
# Sets third coordinate of a vertex for top and bottom interfaces
def set_vert_coord_bdry(ugrid, interface) :
    # Get point coordinate array
    coords = vtk_to_numpy(ugrid.GetPoints().GetData())

    if interface == 0 :    # ground
        coords[:,2] = 1.0 
    elif interface == nz : # top of atmosphere
        coords[:,2] = 0.0 

    # Update vertical coordinate for this layer
    ugrid.GetPoints().SetData(numpy_to_vtk(coords))
    

#########################################################################################################################################
# Functions to convert between images and numpy arrays
def average_vti_images(vti_images) :
    """
    Average the data arrays of a list of VTK image data objects.
    
    Args:
        vti_images (list): A list of VTK image data objects to average.
    
    Returns:
        vtk.vtkImageData: A new VTK image data object containing the averaged data arrays.
    """
    # Get the number of images
    num_images = len(vti_images)
    
    # Get the first image's point data
    first_image = vti_images[0]
    point_data = first_image.GetPointData()

    # Create an empty list to hold the accumulated data arrays
    accumulated_arrays = {}
    
    # Iterate over the data arrays in the first image to initialize the accumulators
    for i in range(point_data.GetNumberOfArrays()):
        array_name = point_data.GetArrayName(i)
        array = point_data.GetArray(i)
        
        # Initialize the accumulator with zeros
        # Make sure the accumulator has the correct shape (num_tuples, num_components)
        num_tuples = array.GetNumberOfTuples()
        num_components = array.GetNumberOfComponents()
        accumulated_arrays[array_name] = np.zeros((num_tuples, num_components))

    # Accumulate the data from each image
    for image in vti_images:
        point_data = image.GetPointData()
        
        for i in range(point_data.GetNumberOfArrays()):
            array_name = point_data.GetArrayName(i)
            array = point_data.GetArray(i)
            
            # Accumulate the array values tuple by tuple
            for j in range(array.GetNumberOfTuples()):
                accumulated_arrays[array_name][j] += np.array(array.GetTuple(j))

    # Create a new vtkImageData to hold the averaged arrays
    averaged_image = vtk.vtkImageData()
    averaged_image.SetDimensions(first_image.GetDimensions())
    averaged_image.SetSpacing(first_image.GetSpacing())
    averaged_image.SetOrigin(first_image.GetOrigin())
    
    # Add the averaged arrays back to the vtkImageData
    for array_name, accumulated_array in accumulated_arrays.items():
        # Average the accumulated data
        averaged_array = accumulated_array / num_images
        
        # Create a VTK array and fill it with the averaged data
        vtk_array = vtk.vtkFloatArray()
        vtk_array.SetName(array_name)
        vtk_array.SetNumberOfTuples(averaged_array.shape[0])
        vtk_array.SetNumberOfComponents(averaged_array.shape[1])
        
        for idx in range(averaged_array.shape[0]):
            vtk_array.SetTuple(idx, tuple(averaged_array[idx]))
        
        averaged_image.GetPointData().AddArray(vtk_array)

    return averaged_image


def read_vti_images(file_list) :
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


def time_mean() :
    # Computes time means of meridional and zonal averages and vertical profile
    zonal = avg_images("zonal")
    merid = avg_images("merid")

    dims       = zonal.GetDimensions() 
    point_data = zonal.GetPointData()
    num_attrs  = point_data.GetNumberOfArrays()

    data_names       = []
    vertical_profile = []
    for i in range(num_attrs) :
        data_names.append(point_data.GetArrayName(i))
        array = vtk_to_numpy(point_data.GetArray(i))
        array = array.reshape((vert_dim, lat_dim))
        profile = np.mean(array,axis=(1)).tolist()
        vertical_profile.append(profile)

    vertical_profile = np.array(vertical_profile).T
    
    with open(run+"_profile.csv", 'w', newline='') as file :
        writer = csv.writer(file)
        writer.writerow(data_names)
        for row in vertical_profile :
            writer.writerow(row)

        
def avg_images(files) :
    # Averages all image files containing files in name (e.g. zonal or merid)
    
    file_list = []
    for t in range (t1, t2+1) :
        file_list.append(run+'_'+str(t).zfill(4)+"_"+files+".vti")
    
    vti_images = read_vti_images(file_list)
    avg = average_vti_images(vti_images)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(run+"_"+files+"_mean.vti")
    writer.SetInputData(avg)
    writer.Write()

    return avg
    

def transform_to_lonlat(t) :
    # Transform spherical data to lonlat vtk data

    script_name = "xyz2lonlat.py"
    arg1 = run
    arg2 = str(Jmin)
    arg3 = str(Jmax)
    arg4 = '1'
    arg5 = str(nz)
    arg6 = str(t)
    arg7 = str(t)
    arg8 = 'y' # use Delaunay2D filter to remove gaps
    subprocess.run(['python3', script_name, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8])

    
def mean (data) :
    # Computes mean of data

    mean = n = 0
    for x, in zip(data):
        n    += 1
        mean += (x - mean) / n
        
    return mean

def covariance(data1, data2) :
    # Stable one-pass covariance of data1, data2
    # returns sample covariance, and means of both variables
    
    meanx = meany = C = n = 0
    for x, y in zip(data1, data2):
        n += 1
        dx     = x - meanx
        meanx += dx / n
        meany += (y - meany) / n
        C     += dx * (y - meany)

    population_covar = C / n
    
    # Bessel's correction for sample variance
    sample_covar = C / (n - 1)

    return sample_covar, meanx, meany


def merge_mean(mean1, n1, mean2, n2):
    """
    Merge a single mean from two datasets with different means.

    Parameters:
        n1, n2 Sample sizes of the two datasets.
        mean1  Mean in dataset 1.
        mean2  Mean in dataset 2.

    Returns:
        Merged mean values.

    """
    n = n1 + n2
    
    delta = mean1 - mean2

    merged_mean = (n1 * mean1 + n2 * mean2) / n
    
    return merged_mean

def merge_covariance(cov1, mean1_x, mean1_y, n1,
                     cov2, mean2_x, mean2_y, n2):
    """
    Merge a single covariance from two datasets with different means.

    Parameters:
        n1, n2           Sample sizes of the two datasets.
        mean1_x, mean1_y Means of variables X and Y in dataset 1.
        mean2_x, mean2_y Means of variables X and Y in dataset 2.
        cov1, cov2       Covariance of (X, Y) in datasets 1 and 2.

    Returns:
        Merged covariance and mean values.

    """
    n = n1 + n2
    
    delta_x = mean1_x - mean2_x
    delta_y = mean1_y - mean2_y

    merged_cov = (n1 * cov1 + n2 * cov2 + (n1 * n2) / (n1 + n2) * delta_x * delta_y) / n

    merged_mean_x = (n1 * mean1_x + n2 * mean2_x) / n
    merged_mean_y = (n1 * mean1_y + n2 * mean2_y) / n
    
    return merged_cov, merged_mean_x, merged_mean_y


def compute_covar(pnt_data, var1, var2) :
    # Computes covariance statistics for variables named var1, var2
    # returns covariance and means of each variable

    x = vtk_to_numpy(pnt_data.GetArray(var1))
    y = vtk_to_numpy(pnt_data.GetArray(var2))

    x = x.reshape((vert_dim, lat_dim, lon_dim))
    y = y.reshape((vert_dim, lat_dim, lon_dim))
    
    covar = np.zeros((vert_dim,lat_dim))
    meanx = np.zeros((vert_dim,lat_dim))
    meany = np.zeros((vert_dim,lat_dim))
        
    for i in range(vert_dim):
        for j in range(lat_dim):
            (covar[i,j], meanx[i,j], meany[i,j]) = covariance(x[i,j,:],y[i,j,:])

    return covar, meanx, meany


def compute_mean_density (pnt_data) :
    # Computes mean density
    
    T    = vtk_to_numpy(pnt_data.GetArray("Temperature"))
    Ps   = vtk_to_numpy(pnt_data.GetArray("Ps"))
    P_Ps = vtk_to_numpy(pnt_data.GetArray("P/Ps"))

    T    =    T.reshape((vert_dim, lat_dim, lon_dim))
    Ps   =   Ps.reshape((vert_dim, lat_dim, lon_dim))
    P_Ps = P_Ps.reshape((vert_dim, lat_dim, lon_dim))

    Rho  = P_Ps * Ps / (Rd * T) # density
    
    Rho_mean = np.zeros((vert_dim,lat_dim))
    for i in range(vert_dim):
        for j in range(lat_dim):
            Rho_mean[i,j] = mean(Rho[i,j,:])

    return Rho_mean

def add_scalar_data(data, N, name, img) :
    # Adds scalar data of total size N to vtk img with given name
    
    attr = numpy_to_vtk(data.reshape(N))
    attr.SetName(name)
    
    img.GetPointData().AddArray(attr)
    
#########################################################################################################################################
#    Main program
#########################################################################################################################################

if (len(sys.argv) < 8) :
    print("""
    Use: python lonlat_to_3D.py run compressible Jmin Jmax nz t1 t2 lon_min lon_max lat_min lat_max vert_min vert_max
    
    Generates a 3D data files, zonal/meridional projections and vertical profiles from a series of layers in directory folder.
    
    Required input parameters:
      run          = prefix name of files (run name)
      compressible = y (compressible) or n (incompressible) simulation
      Jmin         = minimum level
      Jmax         = maximum level
      nz           = number of vertical layers
      t1           = first time
      t2           = last time
    
    Optional parameters:
      lon_min      = minimum longitude
      lon_max      = maximum latitude
      lat_min      = minimum longitude
      lat_max      = maximum latitude
      vert_min     = minimum vertical coordinate in (0,1)
      vert_max     = maximum vertical coordinate in (0,1)

    Example: python lonlat_to_3D.py SimpleJ5J7Z30 y 5 7 30 1 365
    
    Saves the following types of data files:
    run_tttt.vtk            3D unstructured  (lon,lat,P/Ps) 3D vtk data
    run_tttt.vti            3D uniform      (lon,lat,P/Ps) 3D image data
    run_tttt_zonal.vti      2D uniform      (lat,P/Ps)     zonally averaged image data
    run_tttt_merid.vti      2D uniform      (lon,P/Ps)     meridionally averaged image data
    run_tttt_zonal_mean.vti 2D uniform      (lat,P/Ps)     zonally averaged image data averaged over times [t1,t2]
    run_tttt_merid_mean.vti 2D uniform      (lon,P/Ps)     meridionally averaged image data averaged over times [t1,t2]
    run_statistics_mean.vti 2D uniform      (lat,P/Ps)     statistics (temperature variance, eddy momentum flux,
                                                           eddy heat flux, eddy kinetic energy)
    run_tttt.csv            1D                             vertical profiles averaged over the sphere

    3D data has dimensions N x N/2 x nz.  The vertical coordinate is P/Ps.
    """)
    exit(0)
else :
    print("Input parameters = ", sys.argv[1:])

# Input parameters
run          = sys.argv[1]
compressible = sys.argv[2]
Jmin         = int(sys.argv[3])
Jmax         = int(sys.argv[4])
nz           = int(sys.argv[5])
t1           = int(sys.argv[6])
t2           = int(sys.argv[7])

# Dimensions (optional)
if len(sys.argv) > 9 :
    lon_min  = float(sys.argv[8])
    lon_max  = float(sys.argv[9])
    lat_min  = float(sys.argv[10])
    lat_max  = float(sys.argv[11])
    vert_min = float(sys.argv[12])
    vert_max = float(sys.argv[13])

# Grid dimensions
N        = int(np.sqrt(2 * (10*4**Jmax + 2)))
lat_dim  = int(N/2)
lon_dim  = 2*lat_dim
vert_dim = nz

dtheta_min = 360.0/lon_dim
dtheta_max = dtheta_min * 2**(Jmax - Jmin)

lon_min  = -180.0 + dtheta_max
lon_max  =  180.0 - dtheta_max
lat_min  =  -90.0 + dtheta_max
lat_max  =   90.0 - dtheta_max
vert_min =    0.0 
vert_max =    1.0

Ngrid    = lat_dim * lon_dim * vert_dim                                 
Nmerid   = lon_dim * vert_dim
Nzonal   = lat_dim * vert_dim

# Physical constants
Rd          = 287.0  # gas constant
H           = 4000.0 # mean depth
ref_density = 1030.0 # incompressible case
vert_scale  = 1.0

vert_min = vert_scale * vert_min
vert_max = vert_scale * vert_max

# Initialize statistics variables
Ntot      = 0
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

print("\nInterpolating to uniform", lon_dim, "x", lat_dim, "x", vert_dim, "grid")

for t in range (t1, t2+1):
    print("    processing time ", t)
    
    transform_to_lonlat(t) # compute lonlat projections

    vtp_series = []
    for z in range (1, nz+1):
        vtp_series.append(run+'_tri_lonlat_'+str(z).zfill(3)+'_'+str(t).zfill(4)+".vtp")
            
    cell3d = Cell3D(vtp_series)
    cell3d.construct_3Dimage()

# Compute mean over all times
time_mean()

#########################################################################################################################################
