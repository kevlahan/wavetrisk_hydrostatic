#!/usr/bin/env python3

import sys
import os
import math
import vtk


def read_vtp(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def write_vti(image, filename):
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(image)
    writer.Write()


def make_cell_center_sampling_grid(polydata):
    xmin, xmax, ymin, ymax, zmin, zmax = polydata.GetBounds()

    Lx = xmax - xmin
    Ly = ymax - ymin

    if Lx <= 0 or Ly <= 0:
        raise RuntimeError("Degenerate x/y bounds in input .vtp file.")

    n_cells_target = polydata.GetNumberOfCells()

    # Choose nx, ny so nx*ny is approximately equal
    # to the number of cells in the original VTP.
    aspect = Lx / Ly

    nx = max(1, int(round(math.sqrt(n_cells_target * aspect))))
    ny = max(1, int(round(n_cells_target / nx)))

    dx = Lx / nx
    dy = Ly / ny

    # Grid of sampling points located at future VTI cell centres.
    centers = vtk.vtkImageData()
    centers.SetOrigin(xmin + 0.5 * dx,
                      ymin + 0.5 * dy,
                      0.0)
    centers.SetSpacing(dx, dy, 1.0)
    centers.SetDimensions(nx, ny, 1)

    return centers, nx, ny, dx, dy, xmin, ymin


def make_output_cell_grid(nx, ny, dx, dy, xmin, ymin):
    image = vtk.vtkImageData()

    image.SetOrigin(xmin, ymin, 0.0)
    image.SetSpacing(dx, dy, 1.0)

    # Dimensions = number of points
    # => nx*ny cells
    image.SetDimensions(nx + 1, ny + 1, 1)

    return image


def convert_cell_data_to_point_data(polydata):
    """
    vtkProbeFilter samples source point data.
    Convert any cell data to point data first.
    """

    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(polydata)
    c2p.PassCellDataOn()
    c2p.Update()

    return c2p.GetOutput()


def copy_probed_point_data_to_output_cell_data(probed_centers,
                                               output_image):
    pd = probed_centers.GetPointData()
    cd = output_image.GetCellData()

    active_scalar_name = None

    for i in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(i)

        if arr is None:
            continue

        name = arr.GetName()

        # Optional: skip validity mask
        if name == "vtkValidPointMask":
            continue

        cd.AddArray(arr)

        if active_scalar_name is None:
            active_scalar_name = name

    if active_scalar_name is not None:
        cd.SetActiveScalars(active_scalar_name)


def interpolate_vtp_to_cell_centered_vti(vtp_filename,
                                         vti_filename):

    polydata = read_vtp(vtp_filename)

    if polydata.GetNumberOfPoints() == 0:
        raise RuntimeError("Input .vtp contains no points.")

    source = convert_cell_data_to_point_data(polydata)

    centers, nx, ny, dx, dy, xmin, ymin = \
        make_cell_center_sampling_grid(source)

    probe = vtk.vtkProbeFilter()
    probe.SetInputData(centers)
    probe.SetSourceData(source)
    probe.Update()

    probed_centers = probe.GetOutput()

    output_image = make_output_cell_grid(nx, ny,
                                         dx, dy,
                                         xmin, ymin)

    copy_probed_point_data_to_output_cell_data(
        probed_centers,
        output_image
    )

    print(f"Input VTP cells : {polydata.GetNumberOfCells()}")
    print(f"Output VTI cells: {nx} x {ny} = {nx*ny}")
    print(f"dx = {dx:g}")
    print(f"dy = {dy:g}")

    write_vti(output_image, vti_filename)

    print(f"Wrote {vti_filename}")

if __name__ == "__main__":

    if len(sys.argv) == 1:
        vtp_files = sorted(
            f for f in os.listdir(".")
            if f.lower().endswith(".vtp")
        )

        if not vtp_files:
            raise RuntimeError("No .vtp files found in current directory.")

        for input_file in vtp_files:
            root, _ = os.path.splitext(input_file)
            output_file = root + ".vti"

            print()
            print(f"Converting {input_file} -> {output_file}")

            interpolate_vtp_to_cell_centered_vti(
                input_file,
                output_file
            )

    elif len(sys.argv) in (2, 3):
        input_file = sys.argv[1]

        root, ext = os.path.splitext(input_file)

        if ext.lower() != ".vtp":
            raise RuntimeError("Input file must be a .vtp file.")

        if len(sys.argv) == 3:
            output_file = sys.argv[2]
        else:
            output_file = root + ".vti"

        interpolate_vtp_to_cell_centered_vti(
            input_file,
            output_file
        )

    else:
        print("Usage:")
        print("  python vtp_to_cell_vti.py")
        print("  python vtp_to_cell_vti.py input.vtp [output.vti]")
        sys.exit(1)
