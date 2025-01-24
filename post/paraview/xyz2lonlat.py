# Converts xyz data in a vtk file to lonlat data and saves in a new file
# (provides a simple longitude-latitude projection of spherical data)
# Usage: python xyz2lonlat.py input_vtk_file output_vtk_file
# Weiguang Guan (SHARCNET) 2024-04-16

import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# Main program
if (len(sys.argv)<7) :
    print("\nUsage: python xyz2lonlat.py run z1 z2 t1 t2 Delaunay\n")
    print("run      = file base name (without tri and .vtk)")
    print("z1       = first z layer")
    print("z2       = last  z layer")
    print("t1       = first time")
    print("t2       = last time")
    print("Delaunay = y/n (interpolate to Delaunay grid to remove gaps)\n")
    print("output is file_lonlat run_tri_zzz_tttt.vtk or or run_tri_zzz_tttt.vtp (Delaunay=y)\n")
    print("Example:")
    print("python3 xyz2lonlat.py HS_J6J7_dl_240km 1 32 0 28 y")
    exit(0)

file_base = sys.argv[1]
z1        = int(sys.argv[2])
z2        = int(sys.argv[3])
t1        = int(sys.argv[4])
t2        = int(sys.argv[5])
Delaunay  = sys.argv[6]

for z in range (z1, z2+1):
    for t in range (t1, t2+1):
        # Load the input vtk file
        infile  = file_base+"_tri_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        outfile = file_base+"_tri_lonlat_"+str(z).zfill(3)+"_"+str(t).zfill(4)
        print("Transforming file "+infile+".vtk")
        vtkreader = vtk.vtkUnstructuredGridReader()
        vtkreader.ReadAllScalarsOn()
        vtkreader.SetFileName(infile+".vtk")
        vtkreader.Update()

        # Get the unstructed grid data
        unstrctGrid = vtkreader.GetOutput()

        # 1 Get point data and convert them from (x,y,z) to (lon, lat, 0)
        points     = unstrctGrid.GetPoints()
        coords     = vtk_to_numpy(points.GetData())
        num_points = coords.shape[0]

        # Compute radius of sphere
        R = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)

        # Conversion from x,y,z to lon, lat, 0
        coords[:,0] = np.degrees(np.arctan2(coords[:,1], coords[:,0])) # longitude
        coords[:,1] = np.degrees(np.arcsin(coords[:,2] / R))           # latitude
        coords[:,2] = 0.0

        points.SetData(numpy_to_vtk(coords))

        # 2 Get cell array
        cells         = unstrctGrid.GetCells()
        cellformation = vtk_to_numpy(cells.GetData())
        num_cells     = unstrctGrid.GetNumberOfCells()

        # 3 Create point to cell mapping and an array of start indices of cells
        point2cells = [ [] for _ in range(num_points) ] # a list of lists
        startIDs    = [0] * num_cells # start index in cell data "cellformation"

        startID = 0
        for cell in range(num_cells) :    # loop through cells
            startIDs[cell] = startID
            size           = cellformation[startID] # number of vertices
            for i in range(size) :
                pid = cellformation[startID+1+i]
                point2cells[pid].append(cell)
                
            startID = startID + 1 + size

        # 4 Loop over cells, for each which we check if it intersects with the -180/180
        # separation line. If yes, we will move their vertices on one side horizontally
        # so that their lon coordinates are -180. Which side's vertices to be moved 
        # depends on the numbers of vertices of that cell on both sides.
        for cell in range(num_cells) :    # loop through cells
            startID = startIDs[cell]
            #print("startID", startID)
            size = cellformation[startID] # number of vertices
            pids_pos = []
            pids_neg = []
            for i in range(size) :
                pid = cellformation[startID+1+i]
                if (coords[pid,0]<-90 and coords[pid,0]>-180) : # on the negative side
                    pids_neg.append(pid)
                elif (coords[pid,0]>90 and coords[pid,0]<180) : # on the positive side
                    pids_pos.append(pid)

            if (len(pids_pos)!=0 and len(pids_neg)!=0) : # the cell intersects with the line
                if (len(pids_pos)>len(pids_neg)) : # move points on the negative side
                    #print("pos>neg", len(pids_pos), len(pids_neg))
                    for pid in pids_neg :
                        coords[pid,0] = -180.0
                else : # move points on the positive side
                    #print("pos<neg", len(pids_pos), len(pids_neg))
                    for pid in pids_pos :
                        coords[pid,0] = -180.0

        # 5 Split points on the separation line
        last_pid = num_points
        for pid in range(num_points) :
            if (coords[pid,0]==-180.0 or coords[pid,0]==180.0) :
                coords[pid,0] = -180.0 
                points.InsertNextPoint(180.0, coords[pid,1], coords[pid,2])
                coords = np.vstack((coords, [[180.0, coords[pid,1], coords[pid,2]]]))
               
                for cid in point2cells[pid] : # loop over all the cells sharing the point pid
                    # check the cell is on the positive or negative side
                    startID = startIDs[cid]
                    size = cellformation[startID] # number of vertices
                    is_pos = True
                    for i in range(size) :
                        p = cellformation[startID+1+i]
                        if (coords[p,0]!=-180.0) :
                            if (coords[p,0]>0) :
                                is_pos = True
                                cellformation[startID+1+i] = last_pid
                            else :
                                is_pos = False
                            break
                    
                    # replace pid with last_pid for those cells on the positive-lon side
                    if (is_pos) :
                        loc = np.where(cellformation[startID+1:startID+1+size] == pid)
                        cells.ReplaceCellPointAtId(cid, int(loc[0][0]), last_pid)

                        startID = startIDs[cid]
                        size = cellformation[startID] # number of vertices
                        for i in range(size) :
                            if (cellformation[startID+1+i]==pid) :
                                cellformation[startID+1+i] = last_pid
                                break

                last_pid += 1

        # 6 Update the point data
        points.SetData(numpy_to_vtk(coords))

        if Delaunay=='y' :
            # Convert to point data
            cellToPoint = vtk.vtkCellDataToPointData()
            cellToPoint.SetInputData(unstrctGrid)
            cellToPoint.Update()
            output = cellToPoint.GetOutput()

            # Apply the vtkDelaunay2D filter
            delaunay = vtk.vtkDelaunay2D()
            delaunay.SetInputData(output)
            delaunay.Update()
            output = delaunay.GetOutput()

            point2cell = vtk.vtkPointDataToCellData()
            point2cell.SetInputData(output)
            point2cell.Update()
            output = point2cell.GetOutput()

            # Write out structured data
            writer = vtk.vtkXMLPolyDataWriter() 
            writer.SetFileName(outfile+".vtp")
            writer.SetInputData(output)
            writer.Write()
        else:
            writer = vtk.vtkUnstructuredGridWriter()
            writer.SetFileTypeToBinary()
            writer.SetFileName(outfile+".vtk")
            writer.SetInputData(unstrctGrid)
            writer.Write()
         

       

                    

