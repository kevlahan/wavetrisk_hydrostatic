# Average all .vtk, .vti, or .vtp files in current directory.
# Automatically averages both point-data and cell-data arrays.
#
# Assumptions:
#   - directory contains only one of: .vtk, .vti, .vtp
#   - all files have identical geometry/topology
#   - all files have identical point-data and cell-data arrays
#   - array shapes/components match across files

import glob
import numpy as np
import vtk
from vtk.util import numpy_support


# ------------------------------------------------------------
# Detect file type
# ------------------------------------------------------------

vtk_files = sorted(glob.glob("*.vtk"))
vti_files = sorted(glob.glob("*.vti"))
vtp_files = sorted(glob.glob("*.vtp"))

types_present = sum(bool(x) for x in (vtk_files, vti_files, vtp_files))

if types_present > 1:
    raise RuntimeError("Found more than one file type among .vtk, .vti, .vtp. Keep only one type.")

if vti_files:
    files = vti_files
    Reader = vtk.vtkXMLImageDataReader
    Writer = vtk.vtkXMLImageDataWriter
    output_name = "average.vti"
    filetype = "vti"

elif vtp_files:
    files = vtp_files
    Reader = vtk.vtkXMLPolyDataReader
    Writer = vtk.vtkXMLPolyDataWriter
    output_name = "average.vtp"
    filetype = "vtp"

elif vtk_files:
    files = vtk_files
    Reader = vtk.vtkDataSetReader
    Writer = vtk.vtkDataSetWriter
    output_name = "average.vtk"
    filetype = "vtk"

else:
    raise FileNotFoundError("No .vtk, .vti, or .vtp files found.")


def read_dataset(filename):
    reader = Reader()
    reader.SetFileName(filename)

    # Important for legacy .vtk files: read all arrays, not only active/default ones.
    if filetype == "vtk":
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.ReadAllTensorsOn()
        reader.ReadAllNormalsOn()
        reader.ReadAllColorScalarsOn()
        reader.ReadAllFieldsOn()

    reader.Update()
    return reader.GetOutput()


def dataset_signature(data):
    if filetype == "vti":
        return (
            data.GetClassName(),
            tuple(data.GetExtent()),
            tuple(data.GetOrigin()),
            tuple(data.GetSpacing()),
            data.GetNumberOfPoints(),
            data.GetNumberOfCells(),
        )

    return (
        data.GetClassName(),
        data.GetNumberOfPoints(),
        data.GetNumberOfCells(),
    )


def get_array_names(data_attributes):
    return [
        data_attributes.GetArrayName(i)
        for i in range(data_attributes.GetNumberOfArrays())
    ]


def initialize_sums(data_attributes):
    sums = {}

    for name in get_array_names(data_attributes):
        arr_vtk = data_attributes.GetArray(name)
        arr_np = numpy_support.vtk_to_numpy(arr_vtk)
        sums[name] = np.zeros_like(arr_np, dtype=np.float64)

    return sums


def accumulate(data_attributes, sums, filename, association_name):
    current_names = get_array_names(data_attributes)

    if set(current_names) != set(sums.keys()):
        raise RuntimeError(
            f"{filename}: {association_name}-data arrays do not match.\n"
            f"Expected: {sorted(sums.keys())}\n"
            f"Found:    {sorted(current_names)}"
        )

    for name in sums:
        arr_vtk = data_attributes.GetArray(name)

        if arr_vtk is None:
            raise RuntimeError(
                f"{filename}: missing {association_name}-data array '{name}'"
            )

        arr_np = numpy_support.vtk_to_numpy(arr_vtk)

        if arr_np.shape != sums[name].shape:
            raise RuntimeError(
                f"{filename}: {association_name}-data array '{name}' "
                f"has shape {arr_np.shape}, expected {sums[name].shape}"
            )

        sums[name] += arr_np


def replace_with_averages(data_attributes, sums, nfiles):
    # Remove original arrays to avoid duplicate names.
    for name in list(sums.keys()):
        data_attributes.RemoveArray(name)

    # Add averaged arrays.
    for name, total in sums.items():
        avg_np = total / nfiles
        avg_vtk = numpy_support.numpy_to_vtk(avg_np, deep=True)
        avg_vtk.SetName(name)
        data_attributes.AddArray(avg_vtk)

    if sums:
        data_attributes.SetActiveScalars(next(iter(sums.keys())))


# ------------------------------------------------------------
# Read first file and initialize
# ------------------------------------------------------------

first = read_dataset(files[0])

sig0 = dataset_signature(first)
npoints = first.GetNumberOfPoints()
ncells = first.GetNumberOfCells()

point_sums = initialize_sums(first.GetPointData())
cell_sums  = initialize_sums(first.GetCellData())

print(f"Detected file type: .{filetype}")
print(f"Found {len(files)} files")
print(f"Number of points: {npoints}")
print(f"Number of cells:  {ncells}")
print(f"Point-data arrays: {list(point_sums.keys())}")
print(f"Cell-data arrays:  {list(cell_sums.keys())}")


# ------------------------------------------------------------
# Accumulate over files
# ------------------------------------------------------------

for ifile, filename in enumerate(files, start=1):
    print(f"[{ifile}/{len(files)}] {filename}")

    data = read_dataset(filename)

    sig = dataset_signature(data)
    if sig != sig0:
        raise RuntimeError(
            f"{filename}: geometry/topology signature differs.\n"
            f"Expected: {sig0}\n"
            f"Found:    {sig}"
        )

    accumulate(data.GetPointData(), point_sums, filename, "point")
    accumulate(data.GetCellData(),  cell_sums,  filename, "cell")


# ------------------------------------------------------------
# Create output from first file geometry/topology
# ------------------------------------------------------------

output = first.NewInstance()
output.DeepCopy(first)

replace_with_averages(output.GetPointData(), point_sums, len(files))
replace_with_averages(output.GetCellData(),  cell_sums,  len(files))


# ------------------------------------------------------------
# Write averaged file
# ------------------------------------------------------------

writer = Writer()
writer.SetFileName(output_name)
writer.SetInputData(output)

if filetype == "vtk":
    writer.SetFileTypeToBinary()
    # For debugging legacy VTK output, you can instead use:
    # writer.SetFileTypeToASCII()

if writer.Write() != 1:
    raise RuntimeError(f"Failed to write {output_name}")

print()
print(f"✅ Averaged {len(files)} files")
print(f"✅ Output written to '{output_name}'")
