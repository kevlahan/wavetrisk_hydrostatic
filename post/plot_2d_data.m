% Plot 2d data from export_2d
clear all; close all;

machine   = 'mac'
itime     = '030'
itype     = 'vort'
lon_lat   = 1; % Plot longitude - latitude data
zonal_avg = 0; % Plot zonally averaged data
shift     = 1; % shift left boundary to zero longitude
smooth    = 1; % smooth data over two points in each direction

% Extract files
file_base = 'fort.3';
if (strcmp(machine,'if'))
    pathid = '/net/if/1/home/kevlahan/data/jobs/hydrostatic/';
elseif (strcmp(machine,'mac'))
    pathid = '/Users/kevlahan/hydro/';
end
file_tar = ['tar ' 'xvf ' pathid file_base itime '.tgz'];
system(file_tar);

% Load coordinates
lon = load([file_base itime '20']);
lat = load([file_base itime '21']);
P_z = load([file_base itime '22']); % Pressure-based vertical coordinates

if (strcmp(itype,'temp')) % Plot temperature
    %c_scale = 270:3:303; % DCMIP2008c5
    c_scale = 220:10:320; % DCMIP2012c4
    v_title = 'Temperature (K)';
    if (lon_lat)
        s_ll = load([file_base itime '02']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '12']);
    end
elseif (strcmp(itype,'zonal')) % Plot zonal velocity data
    %c_scale = -15:5:50; % DCMIP2008c5
    c_scale = 0:2:20; % DCMIP2012c4
    v_title = 'Zonal velocity (m/s)';
    if (lon_lat)
        s_ll = load([file_base itime '03']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '13']);
    end
elseif (strcmp(itype,'merid')) % Plot meridional velocity data
    %c_scale = -35:5:20;% DCMIP2008c5
    c_scale = -5:1:5;% DCMIP2012c4
    v_title = 'Meridional velocity (m/s)';
    if (lon_lat)
        s_ll = load([file_base itime '04']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '14']);
    end
elseif (strcmp(itype,'geopot')) % Plot geopotential data
    %c_scale = 2400:100:3400; % DCMIP2008c5
    c_scale = 750:100:1500; % DCMIP2012c4
    v_title = 'Geopotential (m)';
    s_ll = load([file_base itime '05']);
    zonal_avg = 0;
elseif (strcmp(itype,'vort')) % Plot relative vorticity data
    s_ll = load([file_base itime '06']);
    c_scale = linspace(-3e-5,3e-5,10);
    v_title = 'Relative vorticity';
    zonal_avg = 0;
end
if (lon_lat)
    if (smooth)
        s_ll = smooth2a(s_ll,2,2);
    end
    plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, 1)
end
if (zonal_avg)
    if (smooth)
        s_zo = smooth2a(s_zo,2,2);
    end
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, 0)
end

fprintf('Minimum value of variable %s = %8.4e\n',itype, min(min(s_ll)));
fprintf('Maximum value of variable %s = %8.4e\n',itype, max(max(s_ll)));

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);

