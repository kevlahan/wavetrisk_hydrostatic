% Plot 2d data from export_2d or log data
clc; clear all;

% Data file details
test_case = 'Held_Suarez'; run_id = 'HS_J6_dl_240km'; run_dir = ''; 
NCAR_topo = false;               % load topography file

%file_base = [run_id '.4.']; % single time
file_base = [run_id '.6.'];  % average

nia       = 0;               % transfer data from niagara if true

Jmin = 6; Jmax = 6;          % minimum and maximum scales (for compression ratio)
Nz = 32;                     % number of vertical layers  (for vertical slice)

% 2d projection options
% 'temp' 'zonal' 'merid' 'geopot' 'vort' 'surf_press' 'ke' 'omega' 
% 'temp_var' 'eddy_mom' 'zonal_var' 'merid_var' 'eddy_ke' 'eddy_heat_flux'

itype     = 'omega';  % field to plot
lon_lat   = true;     % plot longitude - latitude data
zonal_avg = false;    % plot zonally averaged data
shift     = false;     % shift left boundary to zero longitude
smooth    = false;    % smooth data over two points in each direction
lines     = false;    % plot contour lines

if (nia)
    unix ('scp niagara.computecanada.ca:"~/hydro/*c4/DCMIP2012c4_J5J7_var.4.tgz" ~/hydro/DCMIP2012c4/.');
else % local file
   
end

% Log data options:
dt=2; tol_mass=3; tol_temp=4; tol_velo=5; J=6; dof=7; min_mass=8; mass_err=9; balance=10; 
cpu=11; cpudof=12; compression=13; ilog = dof;

% Number of dof on equivalent uniform grid
Nunif = 4 * 10*4^Jmax;

% Total number of degrees of freedom over all scales (this is what is
% counted)
Nmax = 4 * 10*4^Jmin;
for j=Jmin+1:Jmax
    Nmax = Nmax + 4 * 10*4^j;
end

machine = 'mac';
if (strcmp(machine,'if'))
    pathid = ['/net/if/1/home/kevlahan/data/jobs/' test_case '/'];
elseif (strcmp(machine,'mac'))
    pathid = ['/Users/kevlahan/hydro/' test_case '/' run_dir];
end

%set(0,'defaulttextinterpreter','latex')
% Uncompress data files for 2d plots

% Extract files
file_tar = ['tar ' 'xf ' pathid file_base 'tgz'];
disp(['Uncompressing file ' pathid file_base 'tgz']);
system(file_tar);

% Plot 2d data
if (strcmp(itype,'temp_var')||strcmp(itype,'eddy_mom')||strcmp(itype,'eddy_ke')||strcmp(itype,'eddy_heat_flux'))
    lon_lat = 0;
    zonal_avg = 1;
end

% Axis limits
if (strcmp(test_case,'DCMIP2008c5')||strcmp(test_case,'Held_Suarez'))
    if shift
        ax = [0 360 -90 90];
    else
        ax = [-180 180 -90 90];
    end
elseif (strcmp(test_case,'DCMIP2012c4'))
    if (strcmp(itype,'temp')||strcmp(itype,'surf_press'))
        ax = [45 360 0 90];
    else
        ax = [120 270 25 75];
    end
end

% Load coordinates
lon = fread(fopen([file_base '20']),'double');
lat = fread(fopen([file_base '21']),'double');

% Pressure-based vertical coordinates
if (zonal_avg) 
    P_z = fread(fopen([file_base '22']),'double'); 
end

s_ll = 0; s_zo = 0; s_var1 = 0; s_var2 = 0;
if (strcmp(itype,'temp')) % Plot temperature
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 270:3:303;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = 220:10:310;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = 170:10:310;
    end
    v_title = 'Temperature (K)';
    if (lon_lat)
        fid_ll = fopen([file_base '01']);
    end
    if (zonal_avg)
        fid_zo = fopen([file_base '11']);
    end
elseif (strcmp(itype,'temp_var')) % Plot temperature variance
    c_scale = 0:5:50; % Held-Suarez
    v_title = 'Temperature variance (K^2)';
    fid_zo = fopen([file_base '12']);
    zonal_avg = false;
elseif (strcmp(itype,'zonal')) % Plot zonal velocity data
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = -15:5:50;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = -60:5:60;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = -35:5:35;
    end
    v_title = 'Zonal wind (m/s)';
    if (lon_lat)
        fid_ll = fopen([file_base '02']);
    end
    if (zonal_avg)
        fid_zo = fopen([file_base '13']);
    end
elseif (strcmp(itype,'merid')) % Plot meridional velocity data
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = -35:5:20;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = -35:5:35;
    end
    v_title = 'Meridional wind (m/s)';
    if (lon_lat)
        fid_ll = fopen([file_base '03']);
    end
    if (zonal_avg)
        fid_zo = fopen([file_base '14']);
    end
elseif (strcmp(itype,'geopot')) % Plot geopotential data
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 2400:100:3400;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = 750:100:1500;
    end
    v_title = 'Geopotential (m)';
    fid_ll = fopen([file_base '04']);
    zonal_avg = false;
elseif (strcmp(itype,'vort')) % Plot relative vorticity data
    fid_ll = fopen([file_base '05']);
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = -3e-5:1e-5:6e-5;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = -2e-4:5e-5:2e-4;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = linspace(-1e-4,1e-4,100);
    end
    v_title = 'Relative vorticity';
    zonal_avg = false;
elseif (strcmp(itype,'surf_press')) % Plot surface pressure data
    fid_ll = fopen([file_base '06']);
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 800:20:1030;
    elseif (strcmp(test_case,'DCMIP2012c4'))
        c_scale = 930:10:1030;
    end
    v_title = 'Surface pressure (hPa)';
    zonal_avg = false;
elseif (strcmp(itype,'omega')) % Plot vertical velocity
    fid_ll = fopen([file_base '07']);
    c_scale = linspace(-0.1,0.1,100);
     v_title = 'OMEGA at 500 hPa   [Pa/s]';  
elseif (strcmp(itype,'ke')) % Plot kinetic energy
    c_scale = 0:40:480; % Held-Suarez
    v_title = 'Kinetic energy (m^2/s^2)';
    fid_zo = fopen([file_base '15']);
elseif (strcmp(itype,'eddy_mom')) % Plot eddy momentum flux
    c_scale = -80:10:80; % Held-Suarez
    v_title = 'Eddy momentum flux (m^2/s^2)';
    fid_zo = fopen([file_base '16']);
elseif (strcmp(itype,'zonal_var')) % Plot zonal wind variance
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 0:5:100;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = 0:20:300;
    end
    v_title = 'Zonal wind variance(m^2/s^2)';
    fid_zo = fopen([file_base '17']);
elseif (strcmp(itype,'merid_var')) % Plot meridional eddy kinetic energy
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 0:5:100;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = 0:60:600;
    end
    v_title = 'Meridional wind variance (m^2/s^2)';
    fid_zo = fopen([file_base '18']);
elseif (strcmp(itype,'eddy_ke')) % Plot eddy kinetic energy
    if (strcmp(test_case,'DCMIP2008c5'))
        c_scale = 0:5:100;
    elseif (strcmp(test_case,'Held_Suarez'))
        c_scale = 0:40:480;
    end
    v_title = 'Eddy kinetic energy (m^2/s^2)';
    s_zo1 = fopen([file_base '17']) 
    s_zo2 = fopen([file_base '18']);
    fprintf('Mean value of EKE = %8.4e\n', mean(mean(s_zo)));
elseif (strcmp(itype,'eddy_heat_flux')) % Plot zonal eddy heat flux
    c_scale = -22:2:22; % Held-Suarez
    v_title = 'Eddy heat flux (K m/s)';
    fid_zo = fopen([file_base '19']);
end
s_ll = fread(fid_ll,'double'); N = (-3 + sqrt(8*size(s_ll,1) + 1)) / 2;
s_ll = reshape (s_ll, N+1, N/2+1)';

if zonal_avg
    if not(strcmp(itype,'eddy_ke'))
        s_zo = reshape (fread(fid_zo,'double'), N/2, Nz)';
    else
        s_zo = 0.5 * reshape(fread(fid_zo1,'double') + fread(fid_zo1,'double'), N, Nz)';
    end
end

% Plot data
if (lon_lat)
    fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_ll)));
    fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_ll)));
    figure;plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines)
    axis(ax)
    if (NCAR_topo) % extract topography file
        hold on
        topo_file = [run_id '_topo_2D'];
        file_tar = ['tar ' 'xf ' pathid topo_file '.tgz'];
        disp(['Uncompressing topography file ' file_tar]);
        system(file_tar);
        topo = fread(fopen(topo_file),'double'); N = (-3 + sqrt(8*size(topo,1) + 1)) / 2;
        topo = reshape (topo, N+1, N/2+1)';

        topo = circshift(topo,round(size(topo,2)/2),2);
        contour(lon+180,lat,topo, [10 10],'-k','linewidth',2); axis('equal')
    end
end

if (zonal_avg)
    fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_zo)));
    fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_zo)));
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, smooth, lines, 1)
end

% Erase extracted files
file_erase = ['\rm ' file_base '0* ' file_base '2*'];
system(file_erase)
%% Log data plots
% Load log file
%figure;

beg = 1;
%beg = 30000; % Held-Suarez

log_data = load([pathid run_id '_log']);
day = 24;
if ilog == cpudof
    plot(log_data(beg:end,1)/day,log_data(beg:end,cpu)./log_data(beg:end,dof),'k-','linewidth',1.5);
elseif ilog == compression
    plot(log_data(beg:end,1)/day,Nmax./log_data(beg:end,dof),'k-','linewidth',1.5);
else
    plot(log_data(beg:end,1)/day,log_data(beg:end,ilog),'r-','linewidth',1.5);
end

if ilog == dt
    ylab = '$\Delta t$';
elseif ilog == tol_mass
    ylab = 'Mass threshold';
elseif ilog == tol_temp
    ylab = 'Temperature threshold';
elseif ilog == tol_velo
    ylab = 'Velocity threshold';
elseif ilog == J
    ylab = 'J';
elseif ilog == dof
    ylab = '${\cal N}$ (active nodes and edges)';
elseif ilog == min_mass
    ylab = 'Minimum Mass';
elseif ilog == mass_err
    ylab = 'Mass error';
elseif ilog == balance
    ylab = 'Load balance';
elseif ilog == cpu
    ylab = 'CPU time / $\Delta t$ (s)';
elseif ilog == cpudof
    ylab = 'CPU time / ${\cal N}$';
elseif ilog == compression
    ylab = 'Compression ratio';
end

if ilog < cpudof
    axis([0 log_data(end,1)/day 0.9 max(log_data(beg:end,ilog))]);
elseif ilog == cpudof
    axis([0 log_data(end,1)/day 0 max(log_data(beg:end,cpu)./log_data(beg:end,dof))]);
elseif ilog == compression
    axis([0 log_data(end,1)/day 1 1*max(Nmax./log_data(beg:end,dof))]);
end
xlabel('Time (days)');ylabel(ylab);grid on;
set(findall(gcf,'-property','FontSize'),'FontSize',22)