test_case = 'drake'; 
run_id    = '2layer_fill'; 

% Field to plot
set_type = 9;
if set_type == 1
    itype = 'barotropic zonal velocity';
elseif set_type == 2
    itype     = 'barotropic meridional velocity';
elseif set_type == 3
    itype = 'barotropic vorticity';
elseif set_type == 4
    itype = 'layer 1 baroclinic zonal velocity'
elseif set_type == 5
    itype = 'layer 1 baroclinic meridional velocity'
elseif set_type == 6
    itype = 'layer 1 baroclinic vorticity'
elseif set_type == 7
    itype = 'layer 2 baroclinic zonal velocity'
elseif set_type == 8
    itype = 'layer 2 baroclinic meridional velocity'
elseif set_type == 9
    itype = 'layer 2 baroclinic vorticity'
end

lat_min   = -90;
lat_max   =  90;
lon_min   =  0;
lon_max   =  180;

N         = 1024;  % resolution of projection
radius    = 6371.229/6; % radius of Earth

Nx        = N;     % number of points in longitude
Ny        = N/2;   % number of points in latitude
dx        = 2*pi*radius/Nx;
dy        = dx;

smooth    = false;  % smooth data over two points in each direction
shift     = true;   % shift left boundary to zero longitude
lines     = false;  % plot lines

% Load files
run_dir = '';
file_base = [run_id '.4'];
pathid = ['/Users/kevlahan/hydro/' test_case '/' run_dir];
file_tar = ['tar ' 'xf ' pathid file_base '.tgz'];
disp(['Uncompressing file ' pathid file_base '.tgz']);
system(file_tar);

% Load coordinates
lon = load([file_base '.20']); 
lat = load([file_base '.21']);

ax = [lon_min lon_max lat_min lat_max];

if (strcmp(itype,'barotropic zonal velocity'))
    s_ll = load([file_base '.01']);
    c_scale = linspace(-0.8, 1.5, 100);
    v_title = 'Barotropic zonal velocity';
elseif (strcmp(itype,'barotropic meridional velocity'))
    s_ll = load([file_base '.02']);
    c_scale = linspace(-1.1, 1.1, 100);
    v_title = 'Barotropic meridional velocity';
elseif (strcmp(itype,'barotropic vorticity'))
    s_ll = load([file_base '.03']);
    c_scale = linspace(-3.3e-5, 3.3e-5, 100);
    v_title = 'Barotropic vorticity';
elseif (strcmp(itype,'layer 1 baroclinic zonal velocity'))
    s_ll = load([file_base '.04']);
    c_scale = linspace(-0.07, 0.04, 100);
    v_title = 'Baroclinic zonal velocity';
elseif (strcmp(itype,'layer 1 baroclinic meridional velocity'))
    s_ll = load([file_base '.05']);
    c_scale = linspace(-0.06, 0.05, 100);
    v_title = 'Baroclinic meridional velocity';
elseif (strcmp(itype,'layer 1 baroclinic vorticity'))
    s_ll = load([file_base '.06']);
    c_scale = linspace(-2.4e-6, 2.4e-6, 100);
    v_title = 'Layer 1 baroclinic vorticity';
elseif (strcmp(itype,'layer 2 baroclinic zonal velocity'))
    s_ll = load([file_base '.07']);
    c_scale = linspace(-0.1, 0.2, 100);
    v_title = 'Baroclinic zonal velocity';
elseif (strcmp(itype,'layer 2 baroclinic meridional velocity'))
    s_ll = load([file_base '.08']);
    c_scale = linspace(-0.13, 0.16, 100);
    v_title = 'Baroclinic meridional velocity';
elseif (strcmp(itype,'layer 2 baroclinic vorticity'))
    s_ll = load([file_base '.09']);
    c_scale = linspace(-7.2e-6, 7.2e-6, 100);
    v_title = 'Layer 2 baroclinic vorticity';
elseif (strcmp(itype,'free surface'))
    s_ll = load([file_base '.10']);
    c_scale = linspace(-1, 1, 100);
    v_title = 'Free surface';
elseif (strcmp(itype,'internal free surface'))
    s_ll = load([file_base '.11']);
    c_scale = linspace(-350, 1000, 100);
    v_title = 'Internal  freesurface';
elseif (strcmp(itype,'land'))
    s_ll = load([file_base '.12']);
    c_scale = linspace(0, 1, 100);
    v_title = 'Land';
end

fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_ll)));
fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_ll)));
figure
plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines)
axis(ax)
system(['\rm ' run_id '.4*']);