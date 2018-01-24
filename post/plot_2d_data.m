% Plot 2d data from export_2d

clear all; close all;

itime     = '001'
itype     = 'temp'
lon_lat   = 1; % Plot longitude - latitude data
zonal_avg = 0; % Plot zonally averaged data
shift     = 1; % shift left boundary to zero longitude

% Extract files
file_base = 'fort.3';
%pathid = '/net/if/1/home/kevlahan/hydro/';
pathid = '/Users/kevlahan/hydro/';
file_tar = ['tar ' 'xvf ' pathid file_base itime '.tgz'];
system(file_tar);

% Load coordinates
lon = load([file_base itime '20']);
lat = load([file_base itime '21']);
P_z = load([file_base itime '22']); % Pressure-based vertical coordinates

if (strcmp(itype,'temp')) % Plot temperature
    c_scale = 270:3:303;
    v_title = 'Temperature (K)';
    if (lon_lat)
        s_ll = load([file_base itime '02']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '12']);
    end
elseif (strcmp(itype,'zonal')) % Plot zonal velocity data
    c_scale = -15:5:50;
    v_title = 'Zonal velocity (m/s)';
    if (lon_lat)
        s_ll = load([file_base itime '03']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '13']);
    end
elseif (strcmp(itype,'merid')) % Plot meridional velocity data
    c_scale = -35:5:20;
    v_title = 'Meridional velocity (m/s)';
    if (lon_lat)
        s_ll = load([file_base itime '04']);
    end
    if (zonal_avg)
        s_zo = load([file_base itime '14']);
    end
elseif (strcmp(itype,'geopot')) % Plot geopotential data
    c_scale = 2400:100:3400;
    v_title = 'Geopotential (m)';
    s_ll = load([file_base itime '05']);
    zonal_avg = 0; 
end
if (lon_lat)
    plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, 1)
end
if (zonal_avg)
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, 0)
end

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);

