# Averages all .vti files in a directory (they must all have same shape and same data)

import glob
import numpy as np
import vtk
from vtk.util import numpy_support

# ---- Get list of VTI files ----
vti_files = sorted(glob.glob("*.vti"))

if len(vti_files) == 0:
    raise FileNotFoundError("No .vti files found in current directory.")

# ---- Read the first file to get array names and metadata ----
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(vti_files[0])
reader.Update()
image = reader.GetOutput()
point_data = image.GetPointData()

n_arrays = point_data.GetNumberOfArrays()
array_names = [point_data.GetArrayName(i) for i in range(n_arrays)]
npts = image.GetNumberOfPoints()

print(f"Found scalar arrays: {array_names}")

# ---- Initialize dictionary to hold sums for each array ----
sums = {
    name: np.zeros(npts, dtype=np.float32)
    for name in array_names
}
nfiles = 0

# ---- Loop through all files and accumulate data ----
for filename in vti_files:
    reader.SetFileName(filename)
    reader.Update()
    image = reader.GetOutput()
    point_data = image.GetPointData()

    for name in array_names:
        array_vtk = point_data.GetArray(name)
        array_np = numpy_support.vtk_to_numpy(array_vtk)
        sums[name] += array_np

    nfiles += 1

# ---- Average and convert back to VTK arrays ----
output_image = vtk.vtkImageData()
output_image.DeepCopy(image)  # Copy geometry and topology

for name in array_names:
    avg_np = sums[name] / nfiles
    avg_vtk = numpy_support.numpy_to_vtk(avg_np, deep=True)
    avg_vtk.SetName(name)
    output_image.GetPointData().AddArray(avg_vtk)

# ---- Optionally set active scalars (first one) ----
output_image.GetPointData().SetActiveScalars(array_names[0])

# ---- Write the averaged image ----
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName("average.vti")
writer.SetInputData(output_image)
writer.Write()

print(f"✅ Averaged {nfiles} files. Output written to 'average.vti'")
