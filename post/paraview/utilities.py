import numpy as np
import vtk

# Contains various useful functions

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

def rms_int (data, field) :
    # Integrated rms

    data_cells = data.GetCellData()
    
    Calculate (data_cells, field, "np.square(v)", "square")
    
    integrate = vtk.vtkIntegrateAttributes()
    integrate.SetInputData(data)
    integrate.SetDivideAllCellDataByVolume(1)
    integrate.Update()

    integrated = integrate.GetOutput().GetCellData()
    int_square = integrated.GetArray("square").GetTuple(0)[0]

    return(np.sqrt(int_square))
