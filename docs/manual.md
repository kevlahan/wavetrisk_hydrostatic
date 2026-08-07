# WAVETRISK 2.4 
## Users manual and code description
2025-04-17

**Nicholas Kevlahan** 															 

Department of Mathematics and Statistics  
McMaster University							 												 
Hamilton, Canada  
kevlahan@mcmaster.ca	

A dynamically adaptive wavelet-based experimental climate model that provides global three-dimensional atmosphere and ocean models.  `WAVETRISK` solves multilayer shallow water equations with a Lagrangian vertical coordinates Boussinesq form using a mimetic method similar to that of `DYNAMICO`. Local horizontal grid resolution is adapted dynamically at each time step using a multilevel wavelet-based method, providing very high local grid refinement that tracks dynamically important structures.

Matthias Aechtner, Thomas Dubos, Gabrielle Ching-Johnson and Peter Lauritzen have also contributed to the development of WAVETRISK.

There are two basic sub-models:
<pre>
<code>
    WAVETRISK-ATMOSPHERE: parameter compressible = .true.    
    WAVETRISK-OCEAN:      parameter compressible = .false.   
</code>
</pre>  

#### References
* Kevlahan, N K-R & Lemarié, F 2022 WAVETRISK 2.1: An adaptive dynamical core for ocean modelling. *Geosci Model Dev* **15**(17), 6521-6539, [doi.org/10.5194/gmd-15-6521-2022](http:doi.org/10.5194/gmd-15-6521-2022).

* Kevlahan, N K-R 2021 Adaptive wavelet methods for Earth systems modelling. *Fluids* **6** [doi.org/10.3390/fluids6070236](http:doi.org/10.3390/fluids6070236).      

* Kevlahan, N K-R & Dubos, T. 2019 WAVETRISK-1.0: an adaptive wavelet hydrostatic dynamical core. *Geosc Model Dev* **12**, 4901–4921. [doi.org/10.5194/gmd-12-4901-2019](http:doi.org/10.5194/gmd-12-4901-2019).     

* Kevlahan, N K-R, Dubos, T Aechtner, M 2015 Adaptive wavelet simulation of global ocean dynamics using a new Brinkman volume penalization. *Geosci Model Dev* **8** 3891-3909.  [doi.org/10.5194/gmd-8-3891-2015](http:doi.org/10.5194/gmd-8-3891-2015).

* Aechtner, M Kevlahan, N K-R & Dubos, T 2015 A conservative adaptive wavelet method for the shallow-water equations on the sphere. *Q J R Meteorol. Soc.* **141**(690), 1712-1726 [doi:10.1002/qj.2473](http:doi.org/10.1002/qj.2473).            

* Dubos, T & Kevlahan, NK-R 2013 A conservative adaptive wavelet method for the shallow-water equations on staggered grids.  *Q J R Meteorol Soc* **139**, 1997-2020 [doi:10.1002/qj.2097](http:doi.org/10.1002/qj.2097).

### Version history
##### WAVETRISK 2.4
- Separate physics module for use with WAVETRISK-ATMOSPHERE
- Simple dry physics option (Hourdin 1993).
- Redefined the Held-Suarez physics as an option in the physics module.
- Subgrid scale orography parameterization ([Lott and Miller 1997](http:doi.org/10.1002/qj.49712353704), [Japanese Meteorological Agency 2019](http://www.jma.go.jp/jma/jma-eng/jma-center/nwp/outline2019-nwp/pdf/outline2019_03.pdf)).
- Multiscale topography based on the NCAR global topography model.
- WAVETRISK now saves all vertical layers directly as .vtk (paraview) format files  on a single adaptive triangular grid (instead of one grid for each resolution level).
- Post-processing is now handled using dedicated python scripts which use the saved .vtk data files (instead of using matlab).
- Layer data can be easily assembled into a full longitude-latitude-P/Ps three-dimensional uniform grid for visualization of data analysis.
##### WAVETRISK 2.3
- Full multigrid elliptic solver on adaptive grid.
- Removal of implicit horizontal diffusion option.
##### WAVETRISK 2.2
- Dynamic load balancing and automatic checkpointing using `charm++/AMPI`. Removed due to instability and bugs.
##### WAVETRISK 2.1
- Spherical harmonics test case and associated matlab m-file for computing power spectra.
- Ocean modelling test cases.
- Vertical diffusion for ocean modelling using a TKE closure similar to that used in `NEMO`.
- Option of implicit time step for horizontal diffusion.
- Option of semi-implicit barotropic-baroclinic time integration for `WAVETRISK-OCEAN`.
- Adaptive multigrid solver, using scheduled relaxation Jacobi (SRJ) iterations.

##### WAVETRISK 2.0
- 3D hydrostatic extension of original 2D spherical code `WAVETRISK 1.0`.

##### WAVETRISK 1.0
- 2D shallow water equations on the sphere. 
- Implemented in `fortran 2003` and parallelized using `mpi`.  
- Hybrid tree-patch data structure where the patch size is can be chosen to optimize performance (`4x4` is typical).
- Includes load balancing at each checkpoint write.

##### WAVETRISK 0.1
- 2D shallow water equations on the plane. Implemented in `matlab`.

## Contents 
[1. Standard test cases](#markdown-header-1-standard-test-cases)
>[1.1 climate](#markdown-header-11-climate)   
>[1.2 DCMIP2008c5](#markdown-header-12-dcmip2008c5)  
>[1.3 DCMIP1012c4](#markdown-header-13-dcmip1012c4)  
>[1.4 drake](#markdown-header-14-drake)  
>[1.5 flat projection data](#markdown-header-15-flat-projection-data)     
>[1.6 jet](#markdown-header-16-jet)  
>[1.7 make NCAR topo](#markdown-header-17-make-ncar-topo)    
>[1.8 seamount](#markdown-header-18-seamount)     
>[1.9 spherical harmonics](#markdown-header-19-spherical-harmonics)    
>[1.10 tke1d](#markdown-header-110-tke1d)     
>[1.11 upwelling](#markdown-header-111-upwelling)
   
[2. Compiling](#markdown-header-2-compiling)      
>[2.1 Load balancing](#markdown-header-21-load-balancing)        
>[2.2 Optimization](#markdown-header-22-optimization)  
>[2.3 Coarsest grid](#markdown-header-23-coarsest-grid)   
>[2.4 Examples](#markdown-header-24-examples)    

[3. Running](#markdown-header-3-running)  
>[3.1 How to prepare working directory](#markdown-header-31-how-to-prepare-working-directory)  
>[3.2 Checkpointing and restarting](#markdown-header-32-checkpointing-and-restarting)    
>[3.3 mpi examples](#markdown-header-33-mpi-examples)  

[4. Post-processing](#markdown-header-4-post-processing)  
>[4.1 Viewing spherical data with paraview](#markdown-header-41-viewing-spherical-data-with-paraview)   
>[4.2 Data analysis of vtk files using python scripts](#markdown-header-42-data-analysis-of-vtk-files-using-python-scripts)   
>[4.3 Spherical harmonics global and local spectral analysis](#markdown-header-43-spherical-harmonics-global-and-local-spectral-analysis) 

[5. Atmosphere model specifics](#markdown-header-5-atmosphere-model-specifics)   
>[5.1 Physics models for the atmosphere](#markdown-header-51-physics-models-for-the-atmosphere)  
>[5.2 NCAR topography](#markdown-header-52-ncar-topography)    
>[5.3 Subgrid scale orography](#markdown-header-53-subgrid-scale-orography)  

[6. Ocean model specifics](#markdown-header-6-ocean-model-specifics)     
>[6.1 Barotropic-baroclinic mode splitting](#markdown-header-61-barotropic-baroclinic-mode-splitting)  
>[6.2 Vertical turbulent diffusion model and surface fluxes](#markdown-header-62-vertical-turbulent-diffusion-model-and-surface-fluxes)     
>[6.3 Penalization of solid boundaries](#markdown-header-63-penalization-of-solid-boundaries)  

[7. Model details](#markdown-header-7-model-details)     
>[7.1 Prognostic variables](#markdown-header-71-prognostic-variables)   
>[7.2 Parameters](#markdown-header-72-parameters)  
>[7.3 Horizontal coordinates](#markdown-header-73-horizontal-coordinates)   
>[7.4 Vertical coordinates](#markdown-header-74-vertical-coordinates)   
>[7.5 Horizontal grid structure](#markdown-header-75-horizontal-grid-structure)  
>[7.6 Basic grid data types](#markdown-header-76-basic-grid-data-types)  
>[7.7 Grid elements](#markdown-header-77-grid-elements)   
>[7.8 Indexing of grid elements and neighbours](#markdown-header-78-indexing-of-grid-elements-and-neighbours)     
>[7.9 Domain boundary cells](#markdown-header-79-domain-boundary-cells)   
>[7.10 Calculations on adapted grid](#markdown-header-710-calculations-on-adapted-grid)     
>[7.11 Diagnostic variables](#markdown-header-711-diagnostic-variables)  
>[7.12 Horizontal diffusion](#markdown-header-712-horizontal-diffusion)

[8. Paraview visualization](#markdown-header-8-paraview-visualization)     
>[ 8.1 Selecting objects to display](#markdown-header-81-selecting-objects-to-display)  
>[ 8.2 Overlaying topography](#markdown-header-82-overlaying-topography)   
>[ 8.3 Display time](#markdown-header-83-display-time)  
>[ 8.4 Display text](#markdown-header-84-display-text)  
>[ 8.5 Edit color bar](#markdown-header-85-edit-color-bar)  
>[ 8.6 Number format](#markdown-header-86-number-format)  
>[ 8.7 Surface plots](#markdown-header-87-surface-plots)  
>[ 8.8 Slices through 3D data](#markdown-header-88-slices-through-3d-data)  
>[ 8.9 Set axes (clip)](#markdown-header-89-set-axes)  
>[ 8.10 Transform data](#markdown-header-810-transform-data)  
>[ 8.11 Save current view point](#markdown-header-811-save-current-view-point)  
>[ 8.12 Python scripting](#markdown-header-812-python-scripting)  
>[ 8.13 Save figure](#markdown-header-813-save-figure)  
>[ 8.14 Save animation](#markdown-header-814-save-animation)  
>[ 8.15 Save state](#markdown-header-815-save-state)  

<a id="markdown-header-1-standard-test-cases"></a>
## 1. Standard test cases
Test cases are defined in separate subdirectories of `test/` (e.g. `test/climate`).  Each test case consists of three files:
1. The main program (e.g. `climate.f90, drake.f90`). This file specifies parameters and sub-models that are likely to be fixed for the given test case.  It also includes the time stepping loop, initializations and specifies output.
2. A file setting various case-specific properties: `test_case_module.f90` (sets initial conditions, defines the vertical grid, defines and reads/prints parameters for this test case).
3. An example input parameter file: file `test_case.in`. This file is test-case specific and defines parameters that are likely to be modified for different runs.
New test cases can be defined by including a new sub-directory in `~/wavetrisk_hydrostatic/test` and adding the three files above.
<a id="markdown-header-11-climate"></a>
### 1.1 climate  
Climate simulation for global atmosphere using Held and Suarez (1994) or Simple Physics (Hourdin 1993) subgrid scale model. Can also use multiscale NCAR topography. See [Section 5](#markdown-header-5-atmosphere-model-specifics)  for details of the sub-models. 
<a id="markdown-header-12-dcmip2008c5"></a>
### 1.2 DCMIP2008c5  
DCMIP 2008 test case 5: mountain-induced Rossby wave
<a id="markdown-header-13-dcmip1012c4"></a>
### 1.3 DCMIP1012c4  
DCMIP 2012 test case 4: baroclinic instability (Jablonowski and Williamson 2006, QJR Meteorol Soc 2006, 132, 2943–2975).
<a id="markdown-header-14-drake"></a>
### 1.4 drake 
Simplified Drake passage test case on a small planet (inspired by [Ferreira, Marshall and Rose (2011)](https://doi.org/10.1175/2010JCLI3580.1). Used to investigate western boundary current generated turbulence ([Kevlahan and Poulin 2022](http:doi.org/10.1175/JPO-D-21-0318.1)). 
<a id="markdown-header-15-flat-projection-data"></a>
### 1.5 flat projection data   
Post-processing of checkpoint data to calculate flat projection to longitude-latitude grid and compute various diagnostics and statistics. Deprecated in favour of python scripts processing vtk data files.
<a id="markdown-header-16-jet"></a>
### 1.6 jet  
Baroclinic jet test case based on the beta plane configuration in Soufflet et al (Ocean Modelling 98, 36-50, 2016). Used in [Kevlahan and Lemarié (2022)](http:doi.org/10.5194/gmd-15-6521-2022).   
<a id="markdown-header-17-make-ncar-topo"></a>
### 1.7 make NCAR topo
Computes and saves NCAR topography data (see [Section 5.2](#markdown-header-52-ncar-topography)). 
<a id="markdown-header-18-seamount"></a>
### 1.8 seamount  
Seamount test case evaluates pressure gradient error (Beckmann and Haidvogel 1993). An isolated Gaussian bathymetry with flat isopycnals and hybrid sigma coordinates in z (vertical levels not aligned with isopycnals). Used in [Kevlahan and Lemarié (GMD 15 2022)](http:doi.org/10.5194/gmd-15-6521-2022).
<a id="markdown-header-19-spherical-harmonics"></a>
### 1.9 spherical harmonics 
Uses SHTOOLS package to compute global and local spherical harmonics power spectra from checkpoint data interpolated to a uniform grid and projected to a uniform latitude longitude grid ([Wieczorek  and F. J. Simons 2007](http:/doi.org/10.1007/s00041-006-6904-1)). 
<a id="markdown-header-110-tke1d"></a>
### 1.10 tke1d
Implements two 1D test cases for TKE closure scheme for vertical diffusion:
Kato and Phillips (1969) boundary layer thickening due to wind stress forcing (no bottom friction)
Willis and Deardorff (1974) free convection due to surface heat flux (no wind stress, no bottom friction)
see also Zilitinkevich, 1991; Mironov et al., 2000)
<a id="markdown-header-111-upwelling"></a>
### 1.11 upwelling 
Simple wind-driven upwelling/downwelling in a periodic zonal channel. An extension to the sphere of a test case used in ROMS and CROCO. Used in [Kevlahan and Lemarié (GMD 15 2022)](http:doi.org/10.5194/gmd-15-6521-2022).
<a id="markdown-header-2-compiling"></a>
## 2. Compiling
The compile options and their default values are given at the start of the `Makefile`. Note that WAVETRISK requires the gnu compiler collection. WAVETRISK is generally Fortran 2008 compliant. WAVETRISK requires the lapack and netcdf libraries (and shtools/fftw for the spherical harmonics post-processing test case).
<a id="markdown-header-21-load-balancing"></a>
### 2.2 Load balancing   
`MODE=mpi` rebalances the computational load statically at each checkpoint using a simple next fit algorithm.   

<a id="markdown-header-22-optimization"></a>
### 2.3 Optimization  
Default is `OPTIM=2`. Setting `OPTIM=0` adds the `flag -g` and other compiler-dependent checks. `OPTIM=3` or `--ffast-math` is not recommended as it is unstable and causes crashes.
<a id="markdown-header-23-coarsest-grid"></a>
### 2.4 Coarsest grid  
The coarsest grid is set by the compile option `PARAM=param_Jn`, where `n = 4, 5, 6, 7, or 8` is the number of subdivisions of the icosahedron. 

The maximum number of computational cores must be less than or equal to the number of domains, i.e. <code>Ncore <= N_GLO_DOMAIN = 10 2<sup>2 DOMAIN_LEVEL</sup></code>.

Since `DOMAIN_LEVEL = MIN_LEVEL - (PATCH_LEVEL+1)`, larger `MIN_LEVEL` or smaller `PATCH_LEVEL` allows more cores to be used.  For inhomogeneous problems (i.e. unbalanced) it is best to set `PATCH_LEVEL=2` and use a larger `MIN_LEVEL` while for homogeneous (i.e. well-balanced) problems it is more efficient to choose a larger `PATCH_LEVEL`. For example, if `MIN_LEVEL= 7` and `PATCH_LEVEL = 4` then up to 160 cores may be used, while if `PATCH_LEVEL = 2` then up to 2560 cores may be used. These options are set in the files `src/param_Jn.f90`.

The finest allowable grid is set by the parameter max_level in the input file (e.g. `test_case.in`). The finest grid actually used changes dynamically and depends on the solution and the tolerance parameter tol.

There are three options for the coarsest grid, set by the parameter `optimize_grid`:
<pre>
<code>
    optimize_grid = NO_OPTIM    triangle edge bisection
    optimize_grid = XU_GRID     optimized for Laplacian
    optimize_grid = DATA_GRID   optimized for consistency of TRiSK operators
</code>
</pre>
`XU_GRID` (Xu, Int J Comput J, 16(1), 75-93 2006) produces a spherical triangulation that is optimal for the Laplace-Beltrami operator in the sense that the approximation of the discrete mean curvature is exact and the error of the Laplace-Beltrami operator is minimal.

`DATA_GRID` uses an optimized coarse grid with file prefix selected by `grid_type`. The default value `grid_type = HR95JT` uses spherical triangulations (Heikes and Randall 1995, Monthly Weather Rev 123, 1881-1887) that are consistent for the Laplace and flux-divergence operators on the sphere. Laplace, Jacobian and flux-divergence operators have approximately second-order accuracy. Other optimized icosahedral subdivision grids can be added by converting a complete list of vertices to WAVETRISK format using `~/data/vertices2wavetrisk.f90` (run `data/Makefile` to compile). 

The finer grids are then generated by simple edge bisection from the coarsest grid.
<a id="markdown-header-24-examples"></a>
### 2.5 Examples   
Must do `make clean PHYSICS=true` before compiling with the physics sub-model. To compile test case climate using the `Simple Physics` physics model with `J5` as the coarsest grid:
<pre>
<code>
    make TEST_CASE=climate PARAM=param_J5
</code>
</pre>
To compile test case climate using the `Held-Suarez` physics model and NCAR topography:
<pre>
<code>
    make TEST_CASE=climate PARAM=param_J5 TOPO=true
</code>
</pre>
!! Always do "make clean" when compiling a new test case and when you modify the test case !!

<a id="markdown-header-3-running"></a>
## 3. Running  
<a id="markdown-header-31-how-to-prepare-working-directory"></a>
### 3.1 How to prepare working directory  
Initial steps (using test case climate as an example).  

1. cd to the working directory (`mkdir` if necessary), where working directory is where the code will be run.
2. Provide symbolic links to all files required for execution:
<pre>
<code>
    ln -s wavetrisk_hydrostatic/bin/climate
</code>
</pre>

3. Provide a symbolic link to the optimized grids (if `optimize_grid=DATA_GRID`):
<pre>
<code>
    ln -s path-to-grids
</code>
</pre>

The optimized grids are provided in compressed format in `data/grids.xz`.

4. Copy input parameters file and edit as appropriate:
<pre>
<code>
    cp wavetrisk_hydrostatic/test/climate/test_case.in .
</code>
</pre>

If running with bathymetry data (e.g. incompressible tsunami test case) you need to provide a symbolic link to the ETOPO data:
<pre>
<code>
    ln -s ../extern/bathymetry/2arcminutes.smoothed bathymetry (to use 2 arcminute etopo topography data)
</code>
</pre>
See [Section 5.2](#markdown-header-52-ncar-topography) for instruction on how to use NCAR global topography data.
5. Ensure that `gtar` is installed (required for file tar-ing and compression).
<a id="markdown-header-32-checkpointing-and-restarting"></a>
### 3.2 Checkpointing and restarting  
The code checkpoints the data at an interval defined by the test case.  The full state of the simulation is saved by `io.f90/dump_adapt_mpi` and read in again by `io.f90/load_adapt_mpi`.

When a checkpoint is written during a run, the run is restarted immediately from the saved checkpoint. This allows for parallel load balancing and pruning of the adaptive data structure by deleting inactive patches which both occur only during checkpointing. A run can be restarted from a saved checkpoint.

Note that patches that do not include significant cells (i.e. cells in the active zone plus the adjacent zone) are deleted from the data structure ONLY when a checkpoint is saved.  Until then these non-significant cells (labelled with mask `ZERO`) are retained.   A checkpoint+restart therefore changes the adaptive data structure.

Even though non-significant patches are deleted, the adaptive grid always ensures the at least the 4 nearest neighbours of significant cells (i.e. cells outside the adjacent zone) are included in each direction. This allows computation of operators that require at most 5 nearest  neighbours. Any such neighbours that are not significant have the mask value `ZERO`. The inverse wavelet transform onto the adaptive grid ensures these `ZERO` cells have correct data (by interpolation). The trend is computed on the `ADJZONE` mask plus the nearest two neighbouring nodes edges, but only cells in the `ADJZONE` are used for grid adaptation.

Data is written for post processing by routines `io.f90/write_primal` and `io.f90/write_dual` at an interval defined by the test case. Only significant cells are saved (i.e. cells in the active zone plus adjacent zone).
<a id="markdown-header-33-mpi-examples"></a>
### 3.3 mpi examples  
1. With `MODE=mpi`, to run test case climate with input file `climate.in` on a local machine with shared memory 40 cores
<pre>
<code>
    mpirun -n 40 ./climate HS.in
</code>
</pre>
Example slurm script for an mpi job (save in a file, such as up_job.sh):   
<pre>
<code> 
    #!/bin/bash   
    #   
    #SBATCH --account=def-kevlahan        # user account   
    #SBATCH --nodes=16                    # number of nodes: t
    #SBATCH --ntasks-per-node=40 
    #SBATCH --exclusive   
    #SBATCH --time=00-01:00               # maximum running wall clock time (DD-HH:MM)    
    #SBATCH --output=upwelling_J8J10.log  # output file name       
    module load NiaEnv/.2022a gcc openmpi # modules used to compile the executable (up_J8 in this case)   
    srun -n 640 ./up_J8 upwelling.in    # command line to run   
</code>
</pre>
Submit using `sbatch up_job.sh` (where `up_sh.sh` is the name of the slurm script file)   

<a id="markdown-header-4-post-processing"></a>
## 4. Post-processing  
<a id="markdown-header-41-viewing-spherical-data-with-paraview"></a>
### 4.1 Viewing spherical data with paraview  
WAVETRISK saves prognostic variables and some diagnostic variables for each vertical layer on a single triangular cell adaptive grid in *.vtk format. For each data export, vtk files for all layers are compressed into a single `*.vtk.tgz` file. Once uncompressed the vtk data for each layer can be viewed directly in `(x,y,z)` coordinates on paraview as a spherical shell. 

See ~/post/paraview/paraview_notes for tips for visualizing data and making movies using [paraview.org](https://paraview.org).
<a id="markdown-header-42-data-analysis-of-vtk-files-using-python-scripts"></a>
### 4.2 Data analysis of vtk files using python scripts 
In addition to visualization of the spherical data, vtk files are used for data analysis, either using the built-in paraview filters, or by using python code. Some python codes are included in `~/post/paraview`:  
<pre>
<code>
    xyz2lonlat.py.  : Transforms the xyz vtk data to longitude-latitude vtk data without interpolation 
                      (i.e. retaining the triangular cells transformed to lon-lat coordinates).  

                      Usage: python xyz2lonlat.py run z1 z2 t1 t2 Delaunay
                          run      = prefix name of files (without tri)
                          Jmin     = minimum level
                          Jmax     = maximum level
                          z1       = first z layer
                          z2       = last  z layer
                          t1       = first time
                          t2       = last time
                          straight = y/n (straight or jagged longitude edges)
                          Delaunay = y/n (interpolate to Delaunay grid to remove gaps)

                      Output is file_lonlat run_tri_zzz_tttt.vtk or or run_tri_zzz_tttt.vtp (Delaunay=y)
</code>
</pre>

<pre>
<code>
    lonlat_to_3D.py : Generates a 3D data files, zonal/meridional projections and vertical profiles 
                      from a series of layers in directory folder.
 
                      Use: python lonlat_to_3D.py run Jmin Jmax nz t1 t2 lon_min lon_max lat_min lat_max vert_min vert_max
    
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
    
                     Saves the data files:
                          run_tttt.vtk            3D unstructured  (lon,lat,P/Ps) 3D vtk data
                          run_tttt.vti            3D uniform  (lon,lat,P/Ps) 3D image data
                          run_tttt_zonal.vti      2D uniform  (lat,P/Ps)     zonally averaged image data
                          run_tttt_merid.vti      2D uniform  (lon,P/Ps)     meridionally averaged image data
                          run_tttt_zonal_mean.vti 2D uniform  (lat,P/Ps)     zonally averaged image data averaged over [t1,t2]
                          run_tttt_merid_mean.vti 2D uniform  (lon,P/Ps)     meridionally averaged image data averaged over [t1,t2]
                          run_statistics_mean.vti 2D uniform  (lat,P/Ps)     statistics (temperature variance, eddy momentum flux, 
                                                                             eddy heat flux, eddy kinetic energy)
                          run_tttt.csv            1D uniform  (P/Ps)         vertical profiles averaged over the sphere

                     3D data has dimensions N x N/2 x nz.  Vertical coordinate is P/Ps.
</code>
</pre>

<pre>
<code>
    rms.py          : Computes area-integrated rms statistics for various quantities.    

    integrate_ex.py : An example of computing integrated quantities.    

    utilities.py    : Various useful functions.    

    statistics.py   : Taking the difference of two data sets in paraview using the python shell.
</code>   
</pre>
More information about the python codes, including input variables and output, is in included in each script.  

The `matlab` script `post/matlab/plot_vertical_profile.m` plots vertical profiles of diagnostic variables averaged over the sphere using the `.csv` file generated by `lonlat_to_3D.py`.

<a id="markdown-header-43-spherical-harmonics-global-and-local-spectral-analysis"></a>
### 4.3 Spherical harmonics global and local spectral analysis 
The test case spherical_harmonics provides tools for computing spherical harmonics energy spectra. The energy spectra may be computed over the entire sphere, or over a local spherical cap region.  This test case uses the `SHTOOLS` package, which is available at [shtools.github.io/SHTOOLS/index.html](https://shtools.github.io/SHTOOLS/index.html).  The algorithms used in `SHTOOLS` are described in:

* Wieczorek, M. A. and F. J. Simons 2007 Minimum-variance multitaper spectral estimation on the sphere, *J Fourier Anal Appl*, **13**, [doi.org/10.1007/s00041-006-6904-1](http:/doi.org/10.1007/s00041-006-6904-1), 665-692.

* Wieczorek, M. A. and Meschede, M. 2018 SHTools: Tools for Working with Spherical Harmonics. Geochemistry, Geophysics, *Geosystems*, **19**(8), 2574-2592.

`SHTOOLS` requires the libraries `lapack, blas and fftw3`.

An example matlab m-file `spherical_harmonic_analysis.m` is also provided for visualizing the spectra produced by this test case.
 
<a id="markdown-header-5-atmosphere-model-specifics"></a>
## 5. Atmosphere model specifics
The default time integrator is `RK4` and the default adjustment is set to `cfl_num = 0.95`. If `adapt_dt = .true.` the time step is determined based on the smallest value over the grid of `cfl_num * dx / (1/A sum_e(F_e))` w `A` is the hexagon (or pentagon) cell area and `F_e = (|u| + c_s) l_e` is the flux through cell edges, u is the velocity and c_s is the local acoustic wave speed.
<a id="markdown-header-51-physics-models-for-the-atmosphere"></a>
### 5.1 Physics models for the atmosphere
Physics submodels for the atmosphere (e.g. turbulent vertical diffusion of temperature and momentum, seasons, orbit, radiation, convection ) are defined in the directory `~/wavetrisk_hydrostatic/src/physics`.  The atmosphere physics model is turned by including the following line in the test case (see [Section 1.1](#markdown-header-11-climate)):
<pre>
<code>
    physics_model = .true.
</code>
</pre>
which adds an implicit Euler physics time step. There are currently two atmosphere physics models included: `physics/physics_Held_Suarez.f90` (a basic relaxation model defined by Held & Suarez 1995) and `physics/physics_simple.f90` (Hourdin 1992's simple dry physics model). These models are selected using
<pre>
<code>
    physics_type = "Simple"   or
    physics_type = "Held_Suarez"
</code>
</pre>
Each sub-model of Simple Physics can be turned on or off individually using the following flags:
<pre>
<code>
    convecAdj_model  = .true.    convective adjustment module
    diurnal          = .true.    diurnal cycle 
    radiation_model  = .true.    radiation module
    turbulence_model = .true.    vertical diffusion module
</code>
</pre>
All submodels are turned on by default. More physics models could be added to the physics module directory. When using the Simple Physics model the code must be compiled with the flag
<pre>
<code>
    PHYSICS=true
</code>
</pre> 
This flag is automatically set for the `climate` test case.
<a id="markdown-header-52-ncar-topography"></a>
### 5.2 NCAR topography  
Generates smoothed multiscale topography data for WAVETRISK from NCAR topography NetCDF files using `cube_to_target` program that remaps topography data from cubed-sphere grid to target grid (non-adaptive WAVETRISK grid at max_level) grid using rigorous remapping ([Lauritzen, Nair and Ullrich 2015](http:/doi.org/10.5194/gmd-8-3975-2015)).  

The modules `netcdf` and `netcdf-fortran` must be loaded.

You must first generate the multiscale topography for the desired range of levels from `min_level` (as set in the `PARAM` flag) and `max_level` (the maximum resolution level for the test case you want to run). Note that the test case can use any maximum level less than or equal to the `max_level` of the topography.

The complete procedure to generate the multiscale topography is as follows:

1. Compile the code cube_to_target by running `make` in the directory `~/wavetrisk_hydrostatic/topo`. The executable is
<pre>
<code>
    ~/wavetrisk_hydrostatic/bin/cube_to_target  
</code>
</pre>
to generate the `NetCDF` file that provides the surface geopotential `phi_S = z/g` corresponding to the grid data saved in Step 1.  You will need to set the path to the `netcdff` library in `~/wavetrisk_hydrostatic/topo/makefile`, for example
<pre>
<code>
  LDFLAGS = -L/opt/homebrew/Cellar/netcdf-fortran/4.6.2/lib -lnetcdff
</code>
</pre>
The `netcdf` library in `~/wavetrisk_hydrostatic/Makefile` must also be set.

2. Pre-processing of topography data for WAVETRISK test case. Compile the test case `make_NCAR_topo` with `PARAM` set to the coarsest grid resolution (e.g. `PARAM=param_J6`). Then specify the maximum grid resolution in the input file for `make_NCAR_topo` (e.g. `max_level=8`) to generate the WAVETRISK grid coordinates for all levels from `min_level` to `max_level` by sub-sampling. An example input file is:
<pre>
<code>
   run_id          subsample                             ! data set name
   max_level       8                                     ! maximum level of resolution
   topo_data       gmted2010_bedmachine-ncube0540-220518 ! NCAR topoography file 
   smth_scl        30                                    ! smoothing scale for max_leve [km]
   restrict_type   ss                                    ! ss (sub-sampling) or fw (full weighting restriction)
   nsmth_Laplace   0                                     ! number of Laplacian smoothing steps at each level
</code>
</pre>
An example base NCAR topography file for a 3 km resolution cubed sphere is provided in 
<pre>
<code>
   ~/wavetrisk_hydrostatic/data/NCAR_topo/gmted2010_bedmachine-ncube0540-220518.tgz
</code>
</pre>
It is helpful to add symbolic links to the required data files and executables:
<pre>
<code>
    ln -s ~/wavetrisk_hydrostatic/data/gmted2010_bedmachine-ncube0540-220518
    ln -s ~/wavetrisk_hydrostatic/bin/cube_to_target
</code>
</pre>
Note that `make_NCAR_topo` must be run on a single core.

This generates the following topography data file (example for min_level=6, max_level=8 and 30 km smoothing):
<pre>
<code>
    J06J08_030.0km.tgz    topography data for WAVETRISK (including SSO parameters) restricted from max_level to all coarser grids
</code>
</pre>
This multiscale topography data can be used with a test case compiled with `PARAM=param_J6` and `max_level` less than or equal to `8`.

                        
3. The topography_data must have the same `max_level` and domain configuration as the WAVETRISK grid that generated the data in Step 1.  The test case using the NCAR data must be compiled with the flag `TOPO=true` and the modules `netcdff` must be loaded when compiling and running.
<a id="markdown-header-53-subgrid-scale-orography"></a>
### 5.3 Subgrid scale orography model (SSO)
SSO parameterization based on ([Lott and Miller 1997](http:doi.org/10.1002/qj.49712353704), [Japanese Meteorological Agency 2019](http://www.jma.go.jp/jma/jma-eng/jma-center/nwp/outline2019-nwp/pdf/outline2019_03.pdf)). Applicable in the case where there is a large scale separation between the underlying topography data (e.g. NCAR global model) and the finest horizontal grid resolution. This allows the blocking and wave drag effect of unresolved topography to be approximated using the mean and variance of the subgrid scale topography data over a computational grid cell. The `sso` model is activated using the flag
<pre>
<code> 
    sso = .true.
</code>
</pre>
which also requires `NCAR_topo=.true.` and the multiscale NCAR global topography data (See [Section 5.2](#markdown-header-52-ncar-topography)).

A test case using the SSO  must initialize and update at each time step the four statistical quantities mu, theta, gamma and sigma required for the SSO model:
<pre>
<code> 
    sso_param(S_MU)%data(d)%elts
    sso_param(S_THETA)%data(d)%elts
    sso_param(S_GAMMA)%data(d)%elts
    sso_param(S_SIGMA)%data(d)%elts
</code>
</pre>
by calling the subroutine `cal_sso_param` included in the `sso_mod` module. This is done in the routines `/test_case_module/apply_initial_conditions_case` and `/test_case_module/update_case` of the test case (see [Section 1.1](#markdown-header-11-climate) for an example).
<a id="markdown-header-6-ocean-model-specifics"></a>
## 6. Ocean model specifics
<a id="markdown-header-61-barotropic-baroclinic-mode-splitting"></a>
### 6.1 Barotropic-baroclinic mode splitting 
Incompressible (i.e. ocean) simulations can be run with 2D barotropic - 3D baroclinic mode splitting to avoid the stability constraint of the fast external mode with speed c = sqrt(g H). Setting <code>mode_split = .true.</code> chooses mode splitting, which integrates the barotropic component implicitly in time on the slow geostrophic time scale. The computation is stable with CFL_barotropic = sqrt(g H) dt/dx < = 30, and often at larger values. The method is similar to the "implicit free surface" option in MITgcm.  

The splitting uses a semi-implicit theta-method for the barotropic mode, with parameters theta1 (for the external pressure gradient) and theta2 for the barotropic flow divergence.  Theta1 = theta2 = 1 gives a backwards Euler (fully implicit scheme), and theta1 > 0.5, theta2 > 0.5 to avoid the free surface standing wave instability.  The defaults are the NEMO values theta1 = theta2 = 0.55, which ensure stability while avoiding over-damping.

The elliptic equation for the free surface is solved using an adaptive multigrid solver.  Two adaptive multigrid solvers are  provided: Full Multigrid (FMG) using V-cycles and Scheduled Relaxation Jacobi (SRJ) iterations (Adsuara et al J Comput Phys 332, 2017).  In both cases, the  coarsest scales are solved using BiCGSTAB. The default method is FMG,  but SRJ can be selected by setting <code>linear_solver = "SRJ"</code>.

The associated elliptic equation is solved using a simple multigrid scheme on the adaptive grid. The baroclinic and barotropic modes are coupled at each time step. The barotropic part computes an implicit free surface. The barotropic estimate of the location of the free surface is reconciled with the baroclinic estimate using layer dilation (Bleck and Smith 1990, J Geophys Res 95, 3273-3285). Layer dilation means that mass is not exactly conserved in each vertical layer, although it is conserved over all vertical layers. Mode splitting is suitable when the user is not interested in the evolution of the free surface, since it diffuses free surface motion. When the option remap = .false. it is assumed that the density remains constant in each vertical layer, equal to the initial density distribution.  Otherwise, horizontal density gradients are included with Ripa dynamics (e.g. Ripa 1993, Geophys Astrophys Fluid, 70). The drake test case also includes an option that relaxes the layer densities to their initial values. This model is also called IL0 (Beron-Vera 2021, Rev Mex Fis 67), or a "thermal rotating  shallow-water model".

The parameter `level_fill` can be set `> min_level` to make active all `levels <= level_fill`. This sometimes improves the efficiency of the multigrid solver by reducing the size of of the coarsest grid where the elliptic equation is solved exactly.
 
<a id="markdown-header-62-vertical-turbulent-diffusion-model-and-surface-fluxes"></a>
### 6.2 Vertical turbulent diffusion model and surface fluxes    
Eddy diffusion of buoyancy and eddy viscosity for velocity is implemented using a flux-based implicit time integration scheme in the module vert_diffusion_mod (option <code>vert_diffuse = .true.</code>).  Wind stress, bottom friction and buoyancy fluxes are included  via flux (Neumann) boundary conditions. A solar radiation model (as in NEMO) is activated by setting a non-zero value for the solar radiation flux Q_sr. Two options are available: <code>tke_diffuse = .true.</code> uses a turbulent kinetic energy (TKE) closure scheme, similar to that in NEMO ([NEMO book 2015](http:epic.awi.de/id/eprint/39698/1/NEMO_book_v6039.pdf)), <code>tke_diffuse = .false.</code> uses a simple depth-based eddy-viscosity.

<a id="markdown-header-63-penalization-of-solid-boundaries"></a>
### 6.3 Penalization of solid boundaries  
In ocean test cases solid lateral boundaries are modelled using Brinkman volume penalization (see [doi.org/10.5194/gmd-8-3891-2015](http:doi.org/10.5194/gmd-8-3891-2015)).  The permeability parameter eta is set equal to dt (Rasmussen, Cottet & Walther J Comput Phys 230, 2011) and the porosity parameter alpha is set by the user (default value alpha = 0.01). The penalization term is added to the velocity source term in physics_velo_source_case for the test case.

<a id="markdown-header-7-model-details"></a>
## 7. Model details
The numerical approximations are similar to those used in mimetic hydrostatic climate model DYNAMICO [doi.org/10.5194/gmd-8-3131-2015](http:/doi.org/10.5194/gmd-8-3131-2015). The computational grid is based on a hierarchy of bisections of the icosahedron and dual grids of hexagons.
<a id="markdown-header-71-prognostic-variables"></a>
### 7.1 Prognostic variables
The variables are stored in a `Type(Float_Array)` variable `sol` (see below), with components `S_MASS, S_TEMP and S_VELO`.  The code solves the hydrostatic Boussinesq equations in compressible (i.e. atmosphere) or incompressible (i.e. ocean) form.
<pre>
<code>
    S_MASS contains the kinematic/inert pseudo density mu = ref_density dz.

    S_TEMP contains the mass-weighted thermodynamic variable theta, Phi = mu theta.
         Compressible case:           theta = potential temperature.
         Incompressible (ocean) case: theta = buoyancy/g 
         Where buoyancy = - g (density - ref_density)/ref_density = -g drho/g. 
         Total density is therefore ref_density (sol(S_MASS) - sol(S_TEMP)) / sol(S_MASS).

    S_VELO contains the three velocity components U, V, W at edges RT, DG and UP respectively (see below).
</code>
</pre> 
<a id="markdown-header-72-parameters"></a>
### 7.2 Parameters  
The full list of parameters is given in <code>shared.f90</code>, together with their default values.
<a id="markdown-header-73-horizontal-coordinates"></a>
### 7.3 Horizontal coordinates
Cartesian `(x,y,z)` coordinates are used with spherical geometry for lengths and areas on the sphere.  Since the equations are written in vector invariant form, there are no explicit metric terms.   
<a id="markdown-header-74-vertical-coordinates"></a>
### 7.4 Vertical coordinates 
The zero level of the vertical coordinate is at mean sea level and vertical layers are indexed starting at the lowest level in both the atmosphere and ocean models.  The coordinate of the local topography is given by <code>dom%topo%elts(id_i)</code> (or by the test case function surf_geopot in the case of analytical topography), which is positive for orography and negative for bathymetry.  For incompressible (i.e. ocean) simulations the mean depth is defined to be negative, <code>H0 = dom%topo%elts(id_i) < 0</code>, and positive perturbations to the free surface are positive.  The local total depth is <code> H = sum_zlev ( sol(S_MASS,zlev)%data(d)%elts(id_i) ) / ref_density</code> and the free surface perturbation is therefore H + H0.

The vertical coordinates are Lagrangian (i.e. they move as material surfaces) with periodic conservative remapping onto a target (e.g. initial grid).  A collection of different conservative piecewise remapping routines are available in <code>src/remap.f90</code>.
<a id="markdown-header-75-horizontal-grid-structure"></a>
### 7.5 Horizontal grid structure  
The icosahedron is divided into a network of ten regular lozenge grids (with the exception of the two `POLE`). Each lozenge is then divided into <code>N_SUB_DOM = 2<sup>2 DOMAIN_LEVEL</sup> </code> regular sub-domains with the number of subdomains on each processor given by `n_domain(rank+1)`. The lozenges are in turn composed of patches of sub-lozenges, each composed of two triangular (primal) cells.  The lozenge cells each include one hexagon centre `NODE` (for scalars), three `EDGES` (for velocities/fluxes) and 
two triangle centre `NODE`s (for circulation). (See Sections [7.7](#markdown-header-77-grid-elements) and [7.8](#markdown-header-78-indexing-of-grid-elements-and-neighbours).)

There are two special nodes called <code>POLE</code>s. One connects the five lozenge vertices at the top of the network and the other connects the five lozenge grid vertices at the bottom of the network. The network is oriented so the `POLE` nodes are at the geographic north and south poles. The north pole and south poles are stored in the <code>NORTHWEST</code> and <code>SOUTHEAST</code> corner boundaries of two domains. The domains containing poles have <code>grid(d)%pole_master(c/2-2) = .true.</code>, where `c = NORTHWEST</code> or <code>SOUTHEAST`. Applying a routine to the boundaries with options `st = 0, en = 1` applies it to the poles. The boundary choice `(0,1)` should NOT be used for routines that modify variables defined at edges since the poles do not have associated edges (see [Section 7.10](#markdown-header-710-calculations-on-adapted-grid)). The routine apply_to_pole can be used when you do not want to apply a routine only to the poles.

The coarsest level <code>MIN_LEVEL = DOMAIN_LEVEL + PATCH_LEVEL + 1</code> with <code>PATCH_LEVEL>=2</code>. The geometry of this coarsest level may be optimized by reading in a Heikes & Randall (1995) (`optimize_grid = DATA_GRID`) grid of resolution `MIN_LEVEL-1` or using the Xu (2006) smoothing algorithm (<code>optimize_grid = XU_GRID</code>). The total number of domains is 
<pre>
<code>
N_GLO_DOMAIN = N_ICOSAH_LOZENGE N_SUB_DOM = 10 2<sup>2 DOMAIN_LEVEL</sup>
</code>
</pre>
The data structure on each domain `d` is selected as `grid(d)` or `dom=grid(d)`, for example
<pre>
<code>
    grid(d)%ke%elts               selects kinetic energy on domain d
    sol(S_VELO,zlev)%data(d)%elts selects velocity on domain d and vertical layer zlev
</code>
</pre>

The domains are the basic unit of parallelization, and domains are distributed over cores for load balancing. See Sections [2.2](#markdown-header-24-coarsest-grid) and [2.4](#markdown-header-24-coarsest-grid) for details of how the mpi parallelization uses the domain and patch hybrid data structure. 

If <code>optmize_grid = NONE</code> the coarsest grid is just the grid produced by bisecting the icosahedron the `MIN_LEVEL` times. The size of patches is PATCH_LEVEL (e.g. if <code>PATCH_LEVEL = 2</code> then the patches are `4x4`).  The larger <code>PATCH_LEVEL</code> is the fewer levels of tree structure must be 
traversed before reaching a uniform grid (which will in general contain inactive nodes).  Each computation patch element is made up of one node for scalars and three edges for vectors
<code>U, V, W</code> (see also [Section 7.7](#markdown-header-77-grid-elements)).
<a id="markdown-header-76-basic-grid-data-types"></a>
### 7.6 Basic grid data types 
<code>Type(Coord_Array)</code> has components <code>elts(:)%length</code> (where <code>elts</code> is <code>Type(Coord)</code> array of <code>x%y%z</code> coordinates on the sphere for each grid element and length is an integer giving the size of `elts`).

<code>Type(Float_Array)</code> has components <code>elts(:)%length</code> (where <code>elts</code> is a real array of variable values on each grid element on each sub-domain and and length is an integer giving the size of elts). It is used for physical variables such as <code>coriolis, divu, vort</code> etc.

<code>Type(Int_Array)</code> has components <code>elts(:)%length</code> (where <code>elts</code> is an integer array of element indices on each sub-domain and and length is an integer giving the size of <code>elts</code>). It is used for masks, levels and parallel communication variables.

<code>Type(Topo_Array)</code> has components <code>data(:)%bdry_uptodate%pos</code> (where <code>data</code> is a <code>Float_Array, bdry_uptodate</code> is a <code>logical</code> and <code>pos</code> is an
<code>integer</code>) and is used for NCAR topography data.

<code>Type(Float_Field)</code> has components <code>data(:)%bdry_uptodate%pos</code> (where <code>data</code> is a <code>Float_Array, bdry_uptodate</code> is a <code>logical</code> and <code>pos</code> is an
<code>integer</code>) and is used for equation variables such as <code>sol, trend, wav_coeff, sca_coeff</code>.

<code>Type(Domain)</code> has many components defining the grid and is the size of the total number of sub-domains on each processor <code>n_domain(rank+1)</code>.  It is used for the variable <code>grid(:)</code>.

Various other dynamical data types and subroutines for allocating them are defined in <code>dyn_array.f90</code>.   
<a id="markdown-header-77-grid-elements"></a>
### 7.7 Grid elements    
The fundamental data structure element is a lozenge composed of 1 hexagon centre node `H` (stores scalars), 3  edges `RT, DG, UP` (store velocities and fluxes) and 2 triangular nodes `UPLT, LORT` (store circulations):
<pre>
<code>
              ------------ 
             \           / \ 
              \  UPLT   /   \ 
               \       /     \
               UP     DG      \   
                 \   /   LORT  \
                  \ /           \
                   H ------RT---- 
</code>
</pre>

Patch neighbour directions (based on regular coordinates <code>i,j</code>) are labelled as `EAST, NORTHEAST, NORTH, NORTHWEST, WEST, SOUTHWEST, SOUTH, SOUTHEAST`.
Note that within a patch similar notation is used for shifts in the regular grid indices <code>i and j</code> (e.g. shift `i-1, j+1` is denoted `idNW`).
<pre>
<code>
------------- ------------- ------------- 
\              \             \            \
 \              \             \            \
  \              \             \            \
   \  NORTHWEST   \    NORTH    \ NORTHEAST  \  
    \              \             \            \
     \              \             \            \
       -------------  ------------- ------------- 
       \              \             \            \
        \              \             \            \
         \              \             \            \  
          \     WEST     \      0      \    EAST    \ 
           \              \             \            \
            \              \             \            \
              -------------  ------------- ------------- 
              \              \             \            \
               \              \             \            \ 
                \              \             \            \
                 \   SOUTHWEST  \  SOUTH      \  SOUTHEAST \
                  \              \             \            \
                   \              \             \            \
                     -------------   ------------ ------------ 
</code>
</pre>
<a id="markdown-header-78-indexing-of-grid-elements-and-neighbours"></a>
### 7.8 Indexing of grid elements and neighbours  
Quantities (e.g. mass, velocities) are all stored in a single array, whose elements are organized in patches.  Each patch has regular array coordinates <code>(i,j)</code>.

Patch offset array <code>offs(N_BDRY+1)</code> contains the starting index in the single array for the current patch as <code>offs(0)</code>. The starting indices for neighbouring patches are given by <code>offs(NORTH)</code> etc.

Patch dimension array `dims(2, N_BDRY+1)` gives the dimensions of the current patch as `dims(2,0)` and neighbouring patches as `dims(2, NORTH)` etc.

Function `id = idx(i,j,offs,dims)` returns the element index for the hexagon node with coordinates '(i,j)' on the patch selected by 'offs' with dimensions 'dim'.

The components of the grid elements are then found as:
<pre>
<code>
    elts(id+1)        - one hexagon node grid element (scalars such mass, temperature)    
    elts(EDGEid+e+1)  - three edge elements e = RT, DG, UP, (stores vectors such as velocities and fluxes)  
    elts(TRIAGid+t+1) - two triangle grid elements t = LORT, UPLT (quantities defined on triangular grid cells, e.g. circulation)
</code>
</pre>
Wavelet coefficients are stored at the SAME nodes/edges as the `nodes/edges` they characterize. Neighbouring nodes are indexed as follows:
<pre>
<code>
    idE  = idx(i+1,j,  offs,dims) - node to the EAST
    idNE = idx(i+1,j+1,offs,dims) - node to the NORTHEAST
    idN  = idx(i,  j+1,offs,dims) - node to the NORTH

    idW  = idx(i-1,j,  offs,dims) - node to the WEST
    idSW = idx(i-1,j-1,offs,dims) - node to the SOUTHWEST
    idS  = idx(i,  j-1,offs,dims) - node to the SOUTH
</code>
</pre>
Additional variables provide local grid geometry information. For example:
<pre>
<code>
     elts(id+1)%hex_inv.      - inverse of area of hexagon id
     elts(TRIAG*id+t+1)       - areas of triangles t = UPLT, LORT for cell id
     len%elts(EDGE*id+e+1)    - length of triangle edge e = RT, DG, UP of hexagon id
     pedlen%elts(EDGE*id+e+1) - length of hexagon  edge e = RT, DG, UP of hexagon id

     Triangle areas making up hexagon cell id:
     areas%elts(id+1  )%part(1)
     areas%elts(id+1  )%part(2)
     areas%elts(idE+1 )%part(3)
     areas%elts(idNE+1)%part(4)
     areas%elts(idNE+1)%part(5)
     areas%elts(idN+1 )%part(6)
</code>
</pre>
These and other geometric quantities are defined in the module `init_mod`.

Note that scales (different grid resolutions) are referred to as "levels". For adjacent levels, "parents" refers to the coarse level `l` (coarser grid) and "child" (or "children") refers to the fine level `l+1` (finer grid). Adjacent levels always differ by a factor of two in grid resolution.

The subroutine
<pre>
<code>
    call apply_onescale (routine, l, z_null, 0, 1)
</code>
</pre>
applies the subroutine `routine` to a single scale `l`, while the subroutine 
<pre>
<code>
    call apply_interscale (routine, l, z_null, 0, 0)
</code>
</pre>
applies a subroutine `routine` which involves the coarser scale `l` and the finer scale `l-1`.  
<a id="markdown-header-79-domain-boundary-cells"></a>
### 7.9 Domain boundary cells 
To support domain-based parallelization, each domain is extended by `BDRY_THICKNESS` rows and columns of cells associated to the 8 neighbouring domains ("halos" or "ghost cells").

A routine is applied to nodes (scalars) using commands like:
<pre>
<code>
    call apply_bdry (routine, zlev, 0, 1)  (st = 0, en = 1)

</code>
</pre>
where `(0,1)` indicates the routine is applied to the `EAST` and `NORTH` ghost cells in addition to the the actual domain cells. Setting `(0,1)` instead of `(0,0)` automatically applies the routine to the poles since the two poles are included as special nodes in the `en = 1` ghost cells of two domains.

A routine is applied to edges (velocity) using commands like:
<pre>
<code>
    call apply_bdry (routine, zlev, 0, 0)  (st = 0, en = 0)
</code>
</pre>
since all edges are included in the actual domain. Routines modifying edge values should NOT be applied to `(0,1)` since `(0,1)` includes the two poles which do not have associated edges (edge values would therefore be incorrect). Instead, to include poles in routines that modify both nodes and edges (like remap and physics models), the routine apply_no_bdry should be used:
<pre>
<code>
    call apply_no_bdry (routine, zlev)
</code>
</pre>
This routine supplies an integer to the routine that is 1 if it being applied to edges, thus allowing the routine to exclude edge updates.

After applying these routines, the ghost nodes/edges are NOT correct for ghost cells outside `(0,1)` (for scalars) and `(0,0)` (for velocity) for any variables modified by the routine.  The routine should indicate this by setting
<pre>
<code>
    variable%bdry_uptodate = .false.
</code>
</pre>
where variable is the name of the modified variable (e.g. sol). The ghost cells are then corrected by the update routine:
<pre>
<code>
    call update_bdry (variable, l)    to update at level l only
    call update_bdry (variable, NONE) to update at all levels
</code>
</pre>
The update command should be called by all routines requiring correct node values outside (0,1) and correct edge values outside `(0,0)`.

!! The update routine corrects ghost cell values only on `(-(BDRY_THICKNESS-1), BDRY_THICKNESS)`, not `(-BDRY_THICKNESS, BDRY_THICKNESS)`. Therefore, the default value
<pre>
<code>
    BDRY_THICKNESS = 3
</code>
</pre>   
ensures that 2 nearest neighbours are correct for nodes evaluated on `(0,1)` (i.e. ghost cells `-2, 3`), which is sufficient provided no routine requires more than two nearest neighbour nodes/edges in each direction at the same level.
<a id="markdown-header-710-calculations-on-adapted-grid"></a>
### 7.10 Calculations on adapted grid

By default fields are calculated and operators are applied on the `ENTIRE` grid (including at nodes where the result can be obtained by restriction indicated by `mask=12` and adjacent zone nodes indicated by `mask=8` and the results are then over-written by correct values.  (Note that the solution in the adjacent zone at fine scale j+1 are found from values at coarse scale `j` so the values calculated at scale `j+1` are not actually used.)

This means that some operations on the entire grid could produce intermediate overflows, `inf/NaN`, or invalid indices due to incorrect values at these nodes or their neighbours, even though default values have been set to "reasonable" default values when the grid is extended during adaptivity. Functions and subroutines should take this into account.  Similarly, circulation, vorticity and qe (potential vorticity are first computed incorrectly at pentagon points in `step1` (or `cal_vort`) and then corrected in `post_step1` (or `post_vort`).
<a id="markdown-header-711-diagnostic-variables"></a>
### 7.11 Diagnostic variables
Diagnostic variables are derived from the prognostic variables defined in [Section 7.1](#markdown-header-71-prognostic-variables). Diagnostic variables (or quantities) include temperature, zonal velocity, meridional velocity, `OMEGA` (vertical velocity `D_t P = OMEGA [Pa/s]`), vorticity, eddy kinetic energy, eddy momentum flux. Many of these diagnostics are computed in provided routines (e.g. in the data saving routines `write_tri`, or dedicated routines such as `omega_velocity`). Routines for test case specific diagnostic variables are usually provided in the file `test_case_module.f90` for that test case. You can also compute diagnostic quantities from vtk files during post-processing in paraview, python, or matlab (see [Section 4](#markdown-header-4-post-processing)).

For computing quantities integrated over the entire adapted grid, you should use `integrate_tri` (`fun, zlev`), which integrates the function defined by the routine fun over the entire adaptive triangular grid at vertical level `zlev`.  The routine integrate_hex (`fun, zlev, level`) integrates fun on the grid defined at scale level and vertical level `zlev`.  Note that the triangular cells are not overlapping, but the hexagonal cells overlap between levels.  To compute over a portion of the grid, modify fun by including the appropriate mask.
<a id="markdown-header-712-horizontal-diffusion"></a>
### 7.12 Horizontal diffusion
Although WAVETRISK is stable without horizontal diffusion when run non-adaptively, horizontal diffusion of the prognostic variables is required for adaptive runs for both atmosphere and ocean modes.  It is also recommended for non-adaptive runs to remove grid scale vorticity noise.  

The type of diffusion is specified by the parameters `Laplace_sclr` (scalars), `Laplace_divu` (divergent part of velocity) and `Laplace_rotu` (rotational part of velocity) to `0`, `1`, or `2`:
<pre>
<code>
Laplace order
      0          no horizontal diffusion
      1          Laplacian diffusion
      2          bi-Laplacian hyperdiffusion (default for all variables)
</code>
</pre>
The non-dimensional viscosity is specified based on the maximum level of resolution (i.e. finest grid) for each variable using `C_visc(S_MASS)`, `C_visc(S_TEMP)`, `C_visc(S_DIVU)`, `C_visc(S_ROTU)`.  

Note that for stability `C_visc(S_MASS,S_TEMP,S_DIVU) <= (1/9)`<sup>Laplace order</sup> and `C_visc(S_ROTU) <= (1/9/4)`<sup>Laplace order</sup>. 
The default value is `C_visc = 1.5e-3` and `C_visc(S_DIVU) = (1/9)`<sup>Laplace order</sup> (maximum stable for hyperdiffusion) is recommended to guarantee stability when running the climate test cases with a large range of resolution levels.

<a id="markdown-header-8-paraview-visualization"></a>
## 8. Paraview visualization
WAVETRISK saves prognostic variables and some diagnostic variables for each vertical layer in a set of `*.vtk` https://vtk.org format paraview binary files (one for each layer) on initialization and at time intervals specified by `dt_write` (in days). The type of grid is specified by `vtk_grid` which takes the values `"tri"` (triangular grid, default) or `"hex"` (hexagonal grid).  The `tri` format savesa single triangular cell adaptive grid of non-overlapping cells over all active cells at all levels. The `hex` format saves overlapping hexagonal cells for all active cells at all levels. For each data export, vtk files for all layers are compressed into a single `*.vtk.tgz` file. Once uncompressed the vtk data for each layer can be viewed directly in `(x,y,z)` coordinates on paraview as a spherical shell. 
<a id="markdown-header-81-selecting-objects-to-display"></a>
### 8.1 Selecting objects to display
Load `vtk` data files. Click on the window and then in the pipeline browser: toggle `eye` icon to left of object to activate or de-activate its display in window. May need to open "Properties" window.

<a id="markdown-header-82-overlaying-topography"></a>
### 8.2 Overlaying topography
Load `vtk` data files again, click the `eye` icon and select `Mass` variable.  Open color map editor, check "Enable opacity mapping 
for surface", set `Number of table values` to 2 with values 0 and 1. Click the `Gear` icon in `Mapping data` and verify that the opacity for value 0 is 0 and opacity for value 1 is 1.

The `Color transfer function value` for 1 can be set to, e.g., (0.2, 0.2, 0.4) for land areas.

<a id="markdown-header-83-display-time"></a>
### 8.3 Display time
    Filters -> Annotation -> Annotate Time Filter
Need to click on `Apply` to be able to format font and location.
Uses `printf` format (e.g. `%4.2f`). Use shift and scale to convert filenumber to correct units.

<a id="markdown-header-84-display-text"></a>
### 8.4 Display text
    Sources -> Text
Need to click on `Apply` in `Properties` to be able to format font and location.  Font size 22 points is a good choice for text in axes and labels.

<a id="markdown-header-85-edit-color-bar"></a>
### 8.5 Edit color bar
    View -> Color Map Editor
Click on `rainbow` with black e (upper right) to edit color bar legend. Click on `Gear` to get all properties. 
Usually want to unclick `Add Range Labels`. Click on `heart` to change the colour palette: `Rainbow Desaturated` is a good choice in general.

<a id="markdown-header-86-number-format"></a>
### 8.6 Number format
Paraview uses printf formats: e.g. `%1.1f` (standard) or `%1.0e` (scientific notation).
`%#f` always prints decimal point (usually don't want that!).

Use `Interpret values as categories` for simple display of grid resolution levels.

<a id="markdown-header-87-surface-plots"></a>
### 8.7 Surface plots
Ensure you click "apply" in `Properties` to be able to use filters (like contour or smoothing).

    View -> Equalize views to ensure equal window sizes 

    Set axis labels to 20 points, colorbar thickness to 8 and length to 0.6 (in general)

    surface (not slice) 

    View -> Light inspector

In `Key` set `Int` to 1 to avoid dim colors.  Could also adjust `Ele` (elevation) and `Azi` (azimuth) to ensure direct light (e.g. 0, 0 for z view).

<a id="markdown-header-88-slices-through-3d-data"></a>
### 8.8 Slices through 3D data
To view a slice through a 3D dataset (e.g., a horizontal or vertical plane at a specified location):

Load dataset (e.g., `*.vti`), click `Apply` in `Properties`

<pre>
<code>
    Filter -> Common -> Slice
</code>
</pre>
and click open `Eye` in `pipeline`.

In `Properties` panel select `Slice type -> Plane`

In `Plane` parameters panel unclick `Show plane` and select appropriate direction (`X Normal` for zonal slice, `Y Normal` for meridional slice, `Z Normal` for horizontal slice)

Click appropriate `View` direction (`-X` for zonal slice, `-Y` for meridional slice, `-Z` for horizontal slice)
  
Select `depth` in `Origin` boxes (e.g., `[0,0,1` for `P/Ps = 1` for a horizontal slice)

Don't forget to rescale vertical data as appropriate (e.g. by factor 100 in z direction for `P/Ps` data)

Select variable (e.g., `Vorticity`)

In Properties panel check Data Axes Grid, then click Edit and click Gear to get all options. Set font sizes (e.g., 20 for `Title font`, 18 for Axis Labels)
    
Set `Axes` To Label (e.g., `Min-Y`, `Min-Z`)

In `Background` unclick Use `Color Palette` for `Backgroun`d and set background colour to `white` 

Click `Apply`

<a id="markdown-header-89-set-axes"></a>
### 8.9 Set axes
To set bounding axes for a figure you use the `Clip` filter
<pre>
<code>
    Filter -> Common -> Clip
</code>
</pre>
and click on click open `Eye` in `pipeline` to show the `Clip` (unclick other objects in pipeline).  

In `Properties` set `Clip type` to `Box`. 

In `Box parameters` unclick `Show box` and set `Position` (origin) and `Length` appropriately for each direction.  Click `Apply` and select field (e.g. `Vorticity`). 

If the image is blank for a longitude-latitude projection you may need to change the z-direction `Position` and `Length` to ensure the z coordinates of the longitude-latitude projection are included.

<a id="markdown-header-810-transform-data"></a>
### 8.10 Transform data
<pre>
<code>
    Filter -> Alphabetical -> Cell Data to Point Data
</code>
</pre>
To convert `Cell` data (e.g. hexagons) to `Point` data (that can be used by filters like contour).

<a id="markdown-header-811-save-current-view-point"></a>
### 8.11 Save current view point
Click on `camera+` icon to save current view point for future use.

<a id="markdown-header-812-python-scripting"></a>
### 8.12 Python scripting
<pre>
<code>
    View -> Python Shell (opens a python shell to run *.py python scripts)

    Tools -> Start Trace
</code>
</pre>
Can record your modifications as a python script to reuse commands later or in a script.

<a id="markdown-header-813-save-figure"></a>
### 8.13 Save figure
To save the current figure as a `.png` file (i.e. all render views) it is best to work full screen.
<pre>
<code>
    File -> Save Screenshot

    Click Gear to see all options

    Click Save All Views

    Set Image Resolution to actual screen resolution (e.g., 4480 x 2520)

    Set Font Scaling to Scale Fonts Proportionally

    Click Transparent Background

    Set Compression Level to 0 (no compression)

    Save settings and click OK
</code>
</pre>
You may need to adjust font sizes and colour bar positions after checking the `.png` file. Crop the `.png` file to remove extra blanks space using a utility like `Preview`.

<a id="markdown-header-814-save-animation"></a>
### 8.14 Save animation
<pre>
<code>
    File -> Save Animation
</code>
</pre>
Click on `Save All Views` and then save as `.jpeg` (don't click `transparent view`).

Then convert `.jpeg` to `.mov` format usable in keynote (for example):
<pre>
<code>
    ffmpeg -r fps -i base%4d.jpeg -vcodec mpeg4 -q:v 10 base.mov
</code>
</pre>
For example
<pre>
<code>
    ffmpeg -r 12 -i test.%4d.jpeg -vcodec mpeg4 -q:v 10 test.mov
</code>
</pre>
where: `fps` is the number of frames per second, `base` is the base name of `.jpg` files and `%4d` is the frame numbering (e.g. `test.0001.jpeg`, `test.0002.jpeg` etc). 
NB: need to include `-start_number` option if first file does not start at zero 0. 

For example:
<pre>
<code>
    ffmpeg -r 12 -start_number 10 -i test.%4d.jpeg -vcodec mpeg4 --q:v 10 test.mov
</code>
</pre>

<a id="markdown-header-815-save-state"></a>
### 8.15 Save state
To save the complete state of a session (including window layout, dataset names, pipeline, fonts colour bars, etc.) : 
<pre>
<code>
    File -> Save State...
</code>
</pre>
To load a saved state: 
<pre>
<code>
    File -> Load State...
</code>
</pre>