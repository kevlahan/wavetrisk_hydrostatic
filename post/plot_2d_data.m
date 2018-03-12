% Plot 2d data from export_2d
clear all; close all;

machine   = 'if';
t1        = 101; % Start time
t2        = 124; % End time
itype     = 'zonal'; % Options: 'temp' 'zonal' 'merid' 'geopot' 'vort' 'surf_press'
lon_lat   = 0; % Plot longitude - latitude data
zonal_avg = 1; % Plot zonally averaged data
shift     = 1; % shift left boundary to zero longitude
smooth    = 1; % smooth data over two points in each direction
% limits for axis
ax        = [90 200 25 75]; % DCMIP2012c4 vorticity day 7

file_base = 'fort.3';
if (strcmp(machine,'if'))
    pathid = '/net/if/1/home/kevlahan/data/jobs/hydrostatic/';
elseif (strcmp(machine,'mac'))
    pathid = '/Users/kevlahan/hydro/';
end

s_ll = 0; s_zo = 0;
for t = t1:t2
    % Extract files
    itime = num2str(t,'%03i');
    disp(['Loading file ' itime ' of ' num2str(t2,'%03i')]);
    file_tar = ['tar ' 'xf ' pathid file_base itime '.tgz'];
    system(file_tar);
    
    % Load coordinates
    lon = load([file_base itime '20']);
    lat = load([file_base itime '21']);
    P_z = load([file_base itime '22']); % Pressure-based vertical coordinates
    
    if (strcmp(itype,'temp')) % Plot temperature
        %c_scale = 270:3:303; % DCMIP2008c5
        %c_scale = 220:10:320; % DCMIP2012c4
        c_scale = 160:20:300; % Held-Suarez
        v_title = 'Temperature (K)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '02']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '12']);
        end
    elseif (strcmp(itype,'zonal')) % Plot zonal velocity data
        %c_scale = -15:5:50; % DCMIP2008c5
        %c_scale = 0:2:20; % DCMIP2012c4
        c_scale = -30:5:30; % Held-Suarez
        v_title = 'Zonal velocity (m/s)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '03']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '13']);
        end
    elseif (strcmp(itype,'merid')) % Plot meridional velocity data
        %c_scale = -35:5:20;% DCMIP2008c5
        c_scale = -5:1:5;% DCMIP2012c4
        v_title = 'Meridional velocity (m/s)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '04']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '14']);
        end
    elseif (strcmp(itype,'geopot')) % Plot geopotential data
        %c_scale = 2400:100:3400; % DCMIP2008c5
        c_scale = 750:100:1500; % DCMIP2012c4
        v_title = 'Geopotential (m)';
        s_ll = s_ll+load([file_base itime '05']);
    elseif (strcmp(itype,'vort')) % Plot relative vorticity data
        s_ll = s_ll+load([file_base itime '06']);
        %c_scale = -3e-5:1e-5:6e-5;
        c_scale = -1e-5:1e-6:2e-5;
        v_title = 'Relative vorticity';
    elseif (strcmp(itype,'surf_press')) % Plot surface pressure data
        s_ll = s_ll+load([file_base itime '07']);
        c_scale = 800:20:1030; % DCMIP2008c5
        %c_scale = 930:10:1030; % DCMIP2012c4
        v_title = 'Surface pressure';
    end
end

% Average
s_ll = s_ll/(t2-t1+1); s_zo = s_zo/(t2-t1+1);

% Plot data
if (lon_lat)
    fprintf('Minimum value of variable %s = %8.4e\n',itype, min(min(s_ll)));
    fprintf('Maximum value of variable %s = %8.4e\n',itype, max(max(s_ll)));
    if (smooth)
        s_ll = smooth2a(s_ll,2,2);
    end
    plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, 1)
    axis(ax)
end

if (zonal_avg)
    fprintf('Minimum value of variable %s = %8.4e\n',itype, min(min(s_zo)));
    fprintf('Maximum value of variable %s = %8.4e\n',itype, max(max(s_zo)));
    if (smooth)
        s_zo = smooth2a(s_zo,2,2);
    end
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, 0)
end

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);

