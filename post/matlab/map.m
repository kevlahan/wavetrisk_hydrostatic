clear all; clc;
test_case='Held_Suarez';run_id = 'HS_J6'; run_dir = ''; 

pathid = ['/Users/kevlahan/hydro/' test_case '/' run_dir];
% Extract files
topo_file = [run_id '_topo_2D'];
file_tar = ['tar ' 'xf ' pathid topo_file '.tgz'];
disp(['Uncompressing topography file ' topo_file '.tgz']);
system(file_tar);

topo = fread(fopen(topo_file),'double'); N = (-3 + sqrt(8*size(topo,1) + 1)) / 2;
topo = reshape (topo, N+1, N/2+1)'; 

contour(topo, [10 10],'-k','linewidth',2); axis('equal')


if (NCAR_topo) % extract topography file
        file_tar = ['tar ' 'xf ' pathid file_base '.tgz'];
        disp(['Uncompressing file ' file_tar]);
        system(file_tar);

        topo = fread(fopen(topo_file),'double'); N = (-3 + sqrt(8*size(topo,1) + 1)) / 2;
        topo = reshape (topo, N+1, N/2+1)';
    end