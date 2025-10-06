import os
import sys
import fnmatch
import numpy as np
import math
import vtk
import tarfile

# Contains various useful functions

def safe_extract(tar, path=".", members=None):
    # Safely extract files, ensuring no files are extracted outside the target directory.

    for member in tar.getmembers():
        if not os.path.abspath(os.path.join(path, member.name)).startswith(os.path.abspath(path)):
            raise Exception(f"Unsafe extraction attempt: {member.name}")
    tar.extractall(path, members, filter="data")
    
def untar_files (file):
    # Untars time t data

    directory = os.getcwd()
    output_directory = directory
    
    tar_path = os.path.join(directory, file)
    try:
        with tarfile.open(tar_path, 'r:*') as tar:
            safe_extract(tar, path=output_directory)
    except tarfile.TarError as e:
        print(f"    Error extracting {file}: {e}")
    except Exception as e:
        print(f"    Security issue extracting {file}: {e}")

def Calculate (data_cells, field, formula, result_name) :
    # Computes formula applied to the specified field the data_cells
    field_array = data_cells.GetArray(field)

    newDA = vtk.vtkFloatArray()
    newDA.DeepCopy(field_array)
    
    numTuples = newDA.GetNumberOfTuples()
    for i in range(numTuples) :
        v = newDA.GetValue(i)
        newDA.SetValue(i, eval(formula))
        
    newDA.SetName(result_name)
    data_cells.AddArray(newDA)

def rms_sum (arr):
    # Vector rms
    return np.sqrt(np.mean(np.square(arr)))

def delete_files(pattern) :
    # Deletes all files in current directory matching the given pattern.

    directory = os.getcwd()
    
    for filename in os.listdir(directory):
        if fnmatch.fnmatch(filename, pattern):
            file_path = os.path.join(directory, filename)
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")
