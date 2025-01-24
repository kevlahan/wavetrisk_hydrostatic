# Usage: python lonlat_to_3D.py z1 z2 t1 t2 J P_scale
#    
#    Generates a 3D and zonal/meridional projections from a series of vtp layers.
#    
#    Input variables:
#    run     = run
#    z1      = first vertical layer
#    z2      = last  vertical layer
#    t1      = first time
#    t2      = last time
#    J       = scale for the interpolation onto a uniform grid: N/2 x N where N = sqrt(20 4^J)
#    P_scale = scaling factor for P/Ps (for visualization)
#    
#    Saves four data files:
#    run_tttt.vtk            3D unstructured (lon,lat,P) data, where tttt is the time with leading zeros
#    run_tttt.vti            3D uniform      (lon,lat,P) image data
#    run_tttt_zonal_avg.vti  2D uniform      (lat,P) zonally averaged image data
#    run_tttt_merid_avg.vti  2D uniform      (lon,P) merdionally averaged image data
#
#    Data has dimensions N x N/2 x K, where K is the number of vertical layers.
#    The vertical coordinate is kPa.
#
# Elements of *.vtk (also used to compute uniform grid data)
#
#    (1) vtkPoints are a concatenated array of vtkPoints of DS1, DS2, DS3, ..., DNSn.
#
#    (2) vtkCell are constructed by joining the two 2D cells in adjacent layers to produce a prism.
#        (The type of the cell can be vtkWedge, vtkPentagonalPrism, or vtkHexigonalPrism, depending on the type of input 2D cells.)
#
#    (3) vtkCellData are computed as the average of the values from the two 2D cells in adjacent layers.  
#
#
# Author: Weiguang Guan and Nicholas Kevlahan (McMaster University)
# Date  : 2025-01-24
#
# See https://raw.githubusercontent.com/Kitware/vtk-examples/refs/heads/gh-pages/src/Testing/Baseline/Cxx/GeometricObjects/TestLinearCellsDemo.png

import os
from contextlib import suppress
import sys
import numpy as np
import vtk
import csv
import subprocess
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

################################################################################
# This function set the third coordinate of vertex by the average of attribute
# identified by "attr_name" of the cells that share the vertex. For example,
def set_3d_dim(ugrid, attr_name) :
    # Get point coordinate array
    coords     = vtk_to_numpy(ugrid.GetPoints().GetData())
    
    num_points = coords.shape[0]
    weights    = np.zeros(num_points) 

    # Get cell data array with attr_name (which could be 'P/Ps')
    attrs = vtk_to_numpy(ugrid.GetCellData().GetArray(attr_name))
    
    cellformation = vtk_to_numpy(ugrid.GetPolys().GetData())
    num_cells     = ugrid.GetNumberOfCells()

    # Loop through cells 
    startID = 0
    for i in range(num_cells) :
        size    = cellformation[startID]
        pnt_ids = cellformation[(startID+1):(startID+1+size)]
        for pnt in pnt_ids :
            coords[pnt,2] = coords[pnt,2] * weights[pnt] + attrs[i] * P_scale
            weights[pnt] += 1
            coords[pnt,2] /= weights[pnt]
         
        startID = startID + 1 + size

    # Update
    ugrid.GetPoints().SetData(numpy_to_vtk(coords))

    return

################################################################################
# Function to load vtkUnstructuredGrid from a vtk file
def load_dataset(file, attr_to_fill=None) :
    #vtkreader = vtk.vtkUnstructuredGridReader()
    vtkreader = vtk.vtkXMLPolyDataReader() 
    #vtkreader.ReadAllScalarsOn()
    vtkreader.SetFileName(file)
    vtkreader.Update()

    ugrid = vtkreader.GetOutput()
    if (attr_to_fill is not None) :
        set_3d_dim(ugrid, attr_to_fill)

    return ugrid

################################################################################
# Functions to convert between images and numpy arrays
def vtk_image_to_numpy(vtk_image):
    """Convert VTK image data to a NumPy array."""
    dimensions = vtk_image.GetDimensions()  # Get the dimensions of the image
    point_data = vtk_image.GetPointData()
    array = point_data.GetScalars()  # Get the scalar data (image values)
    
    np_array = vtk.util.numpy_support.vtk_to_numpy(array)
    return np_array.reshape(dimensions[2], dimensions[1], dimensions[0])  # Reshape it to 3D

def numpy_to_vtk_image(np_array):
    """Convert a NumPy array to a VTK image."""
    vtk_image = vtk.vtkImageData()
    vtk_image.SetDimensions(np_array.shape[2], np_array.shape[1], np_array.shape[0])
    
    # Flatten the numpy array and set it to the image data
    vtk_array = vtk.util.numpy_support.numpy_to_vtk(np_array.flatten(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_image.GetPointData().SetScalars(vtk_array)
    
    return vtk_image

def average_vtk_images(vtk_images):
    """Average a list of VTK image arrays."""
    numpy_arrays = [vtk_image_to_numpy(img) for img in vtk_images]
    
    # Compute the average of the arrays along the 0th axis (image stacking axis)
    average_array = np.mean(numpy_arrays, axis=0)
    
    # Convert the averaged NumPy array back to VTK
    return numpy_to_vtk_image(average_array)

################################################################################
class Cell3D() :
    def __init__(self, vtp_series):
        global data_names
        self.vtp_series = vtp_series
        self.ugrid = vtk.vtkUnstructuredGrid()  # final result
        
        # Load DS1
        ds1 = load_dataset(vtp_series[0], 'P/Ps')

        points1 = ds1.GetPoints()
        self.ugrid.SetPoints(points1)

        self.num_attrs = ds1.GetCellData().GetNumberOfArrays()

        self.attr_list  = [None] * self.num_attrs     # to store attribute arrays
        self.attr_names = [None] * self.num_attrs

        for i in range(self.num_attrs) :
            attr1 = ds1.GetCellData().GetArray(i)
            self.attr_names[i] = attr1.GetName()

        data_names = self.attr_names

    # Function to construct cells between two layers (ds1 and ds2) and
    # then add them to the unstructured grid (ugrid)
    # Here, ugrid has already had vtkPoints from ds1
    def add_cells(self, ds1, ds2, cell_layer_id=0) :
        # (1) Add vtkPoints
        points = self.ugrid.GetPoints()
        coords = vtk_to_numpy(points.GetData())

        points2 = ds2.GetPoints()
        coords2 = vtk_to_numpy(points2.GetData())

        coords = np.concatenate((coords, coords2))

        points.SetData(numpy_to_vtk(coords))
        self.ugrid.SetPoints(points)

        # (2) Construct vtkCells
        num_points  = points.GetNumberOfPoints()
        num_points2 = ds2.GetNumberOfPoints()

        cellformation = vtk_to_numpy(ds2.GetPolys().GetData())

        num_cells = ds2.GetNumberOfCells()
        startID = 0
        # loop through cells in ds2
        for i in range(num_cells) :
            size     = cellformation[startID]
            pnt_ids  = cellformation[(startID+1):(startID+1+size)]

            pnt_ids1 = pnt_ids + (num_points - 2*num_points2)
            pnt_ids2 = pnt_ids + (num_points -   num_points2)
            
            # 6 points in wedge cell
            pnt_ids  = np.concatenate((pnt_ids1, pnt_ids2))

            # Triangular cells
            self.ugrid.InsertNextCell(vtk.VTK_WEDGE , len(pnt_ids), pnt_ids)
    
            startID = startID + 1 + size

        # (3) Set the cell attributes to be the average of two
        for i in range(self.num_attrs) :
            attr1 = ds1.GetCellData().GetArray(i)
            attr2 = ds2.GetCellData().GetArray(i)
            attr = (vtk_to_numpy(attr1) + vtk_to_numpy(attr2))/2.0

            if (self.attr_list[i] is None) :
                self.attr_list[i] = attr
            else :
                self.attr_list[i] = np.concatenate((self.attr_list[i], attr))

    def construct(self) :
        ds1 = load_dataset(self.vtp_series[0], 'P/Ps')
        print("loading " + self.vtp_series[0])

        for i, file in enumerate(self.vtp_series[1:]) :
            print("loading " + file)
            ds2 = load_dataset(file, 'P/Ps')

            # Construct cells between two layers and add them to ugrid
            self.add_cells(ds1, ds2, (i+1))

            # Last step
            ds1 = ds2

        # Set attributes for the ugrid
        for i in range(self.num_attrs) :
            attr_array = numpy_to_vtk(self.attr_list[i])
            attr_array.SetName(self.attr_names[i])

            self.ugrid.GetCellData().AddArray(attr_array)

        return self.ugrid

    # Returns four elements: unstructured grid, regular grid, and two 2D images as projections
    def construct_3Dimage(self) :
        P_dim = len(self.vtp_series) # max_depth = number of vtk files in the series
        dP = 1.0/P_dim

        P_min   = (0.0 + dP/2) * P_scale
        P_max   = (1.0 - dP/2) * P_scale

        # Construct an unstructured grid
        ugrid = self.construct()

        # Resample ugrid to a regular grid
        ugrid_to_image = vtk.vtkResampleToImage()
        ugrid_to_image.SetInputDataObject(ugrid)
        
        ugrid_to_image.SetUseInputBounds(False)
        ugrid_to_image.SetSamplingDimensions(lon_dim, lat_dim, P_dim)
        ugrid_to_image.SetSamplingBounds(lon_min, lon_max, lat_min, lat_max, P_min, P_max)

        ugrid_to_image.Update()

        rgrid = ugrid_to_image.GetOutput()

        # Generate two projections
        pnt_data = rgrid.GetPointData()

        img1 = vtk.vtkImageData()
        img2 = vtk.vtkImageData()

        vertical_profile = []

        img1.SetDimensions(1, lat_dim, P_dim);
        img2.SetDimensions(lon_dim, 1, P_dim);

        img1.SetSpacing(rgrid.GetSpacing())
        img2.SetSpacing(rgrid.GetSpacing())

        img1.SetOrigin(rgrid.GetOrigin())
        img2.SetOrigin(rgrid.GetOrigin())

        for i in range(self.num_attrs) :
            name = pnt_data.GetArray(i).GetName()

            img3d = vtk_to_numpy(pnt_data.GetArray(i))
            img3d = img3d.reshape((P_dim, lat_dim, lon_dim))

            merid_avg  = np.mean(img3d, axis=2)     # meridional average
            zonal_avg  = np.mean(img3d, axis=1)     # zonal average

            profile = np.mean(img3d, axis=(1,2)).tolist()
            vertical_profile.append(profile) # vertical profile averaged over sphere

            attr1 = numpy_to_vtk(merid_avg.reshape(lat_dim*P_dim))
            attr1.SetName(name)
            attr2 = numpy_to_vtk(zonal_avg.reshape(lon_dim*P_dim))
            attr2.SetName(name)

            img1.GetPointData().AddArray(attr1)
            img2.GetPointData().AddArray(attr2)

        return ugrid, rgrid, img1, img2, vertical_profile
    

################################################################################
# Main program
if (len(sys.argv)<8) :
    print("""
    Use: python lonlat_to_3D.py z1 z2 t1 t2 J P_scale
    
    Generates a 3D vtk data file from a series of layers in directory folder.
    
    Input variables:
    run     = run
    z1      = first vertical layer
    z2      = last  vertical layer
    t1      = first time
    t2      = last time
    t2      = last time
    J       = scale for the interpolation onto a uniform grid: N/2 x N where N = sqrt(20 4^J)
    P_scale = scaling factor for P/Ps (for visualization)
    
    Saves four data files:
    run_tttt.vtk            3D unstructured (lon,lat,P) data, where tttt is the time with leading zeros
    run_tttt.vti            3D uniform      (lon,lat,P) image data
    run_tttt_zonal_avg.vti  2D uniform      (lat,P) zonally averaged image data
    run_tttt_merid_avg.vti  2D uniform      (lon,P) merdionally averaged image data

    Data has dimensions N x N/2 x K, where K is the number of vertical layers.
    The vertical coordinate is kPa.
    """)
    exit(0)

z1      = int(sys.argv[2])
z2      = int(sys.argv[3])
t1      = int(sys.argv[4])
t2      = int(sys.argv[5])
J       = float(sys.argv[6])
P_scale = float(sys.argv[7])

N       = int(np.sqrt(20*4**J))

lat_dim = int(N/2)
lon_dim = 2*lat_dim

dlat = 360.0/lon_dim
dlon = 180.0/lat_dim

lon_min = -180.0 
lon_max =  180.0 - dlon*1.5
lat_min =  -80.0
lat_max =   80.0

# Remove .DS_store to avoid load error
with suppress(OSError):
    os.remove(sys.argv[1]+'/.DS_Store')

#script_path = 'xyz2lonlat.py'
#arguments = ['SimpleJ5J7Z30', '1', '30', '1', '1', 'y']
#print(['python3', script_path] + arguments)
#subprocess.run(['python3', script_path] + arguments, capture_output=True, text=True)

print("\nInterpolating to uniform",lon_dim,"x",lat_dim,"grid\n")

for t in range (t1, t2+1):
    print("\nProcessing time ",t)
    vtp_series = []
    for z in range (z1, z2+1):
        vtp_series.append(sys.argv[1]+'_tri_lonlat_'+str(z).zfill(3)+"_"+str(t).zfill(4)+".vtp")

    cell3d = Cell3D(vtp_series)
    ugrid, rgrid, proj_img1, proj_img2, vertical_profile = cell3d.construct_3Dimage()

    # Save vertical profiles in csv file
    with open(sys.argv[1]+"_"+str(t).zfill(4)+"_profile.csv", 'w', newline='') as file:
        writer = csv.writer(file)
        i = 0
        for row in vertical_profile:
            writer.writerow([data_names[i]]+row)
            i += 1

    # Write out 3D volume
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileTypeToBinary()
    writer.SetFileName(sys.argv[1]+"_"+str(t).zfill(4)+".vtk")
    writer.SetInputData(ugrid)
    writer.Write()

    xmlWriter = vtk.vtkXMLImageDataWriter()

    # Write 3D Cartesian grid data
    xmlWriter.SetFileName(sys.argv[1]+"_"+str(t).zfill(4)+".vti")
    xmlWriter.SetInputData(rgrid)
    xmlWriter.Write()

    # Write zonal projection
    xmlWriter.SetFileName(sys.argv[1]+"_"+str(t).zfill(4)+"_zonal_avg.vti")
    xmlWriter.SetInputData(proj_img1)
    xmlWriter.Write()

    # Write meridional projection
    xmlWriter.SetFileName(sys.argv[1]+"_"+str(t).zfill(4)+"_merid_avg.vti")
    xmlWriter.SetInputData(proj_img2)
    xmlWriter.Write()

