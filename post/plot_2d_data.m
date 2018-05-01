% Plot 2d data from export_2d
clear; close all;

% Test case
test = 'DCMIP2008c5';
%test = 'DCMIP2008c12';
%test = 'Held_Suarez';

machine   = 'mac';
t1        = 2; % Start time
t2        = t1; % End time
% Options: 'temp' 'zonal' 'merid' 'geopot' 'vort' 'surf_press' 'temp_var' 'eddy_mom' 'eddy_ke' 'eddy_heat_flux'
itype     = 'temp';

lon_lat   = 1; % Plot longitude - latitude data
zonal_avg = 0; % Plot zonally averaged data
shift     = 1; % shift left boundary to zero longitude
smooth    = 1; % smooth data over two points in each direction

if (strcmp(itype,'temp_var')||strcmp(itype,'eddy_mom')||strcmp(itype,'eddy_ke')||strcmp(itype,'eddy_heat_flux'))
    lon_lat = 0;
    zonal_avg = 1;
end

% Axis limits
if (strcmp(test,'DCMIP2008c5'))
    ax = [0 360 -90 90];
elseif (strcmp(test,'DCMIP2008c12'))
    ax = [90 200 25 75];
end

N = t2-t1+1; % number of samples

file_base = [test '.3'];
if (strcmp(machine,'if'))
    pathid = '/net/if/1/home/kevlahan/data/jobs/DCMIP2008c5/';
elseif (strcmp(machine,'mac'))
    pathid = '/Users/kevlahan/hydro/DCMIP2008c5/';
end

s_ll = 0; s_zo = 0; s_var1 = 0; s_var2 = 0;
for t = t1:t2
    % Extract files
    itime = num2str(t,'%04i');
    disp(['Loading file ' itime ' of ' num2str(t2,'%04i')]);
    file_tar = ['tar ' 'xf ' pathid file_base itime '.tgz'];
    system(file_tar);
    
    % Load coordinates
    lon = load([file_base itime '20']);
    lat = load([file_base itime '21']);
    P_z = load([file_base itime '22']); % Pressure-based vertical coordinates
    
    if (strcmp(itype,'temp')) % Plot temperature
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = 270:3:303;
        elseif (strcmp(test,'DCMIP2008c12'))
            c_scale = 220:10:320;
        elseif (strcmp(test,'Held_Suarez'))
            c_scale = 160:20:300;
            c_scale2 = 0:0.5:4;
        end
        v_title = 'Temperature (K)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '02']);
        end
        if (zonal_avg)
            s_zo = s_zo + load([file_base itime '12']);
        end
    elseif (strcmp(itype,'temp_var')) % Plot temperature variance
        c_scale = 0:5:40; % Held-Suarez
        v_title = 'Temperature variance (K^2)';
        s_zo = s_zo + load([file_base itime '13']);
    elseif (strcmp(itype,'zonal')) % Plot zonal velocity data
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = -15:5:50; % DCMIP2008c5
        elseif (strcmp(test,'DCMIP2008c12'))
            c_scale = 0:2:20;   
        elseif (strcmp(test,'Held_Suarez'))
            c_scale = -30:5:30;
        end
        v_title = 'Zonal velocity (m/s)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '03']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '14']);
        end
    elseif (strcmp(itype,'merid')) % Plot meridional velocity data
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = -35:5:20;
        elseif (strcmp(test,'DCMIP2008c12'))
            c_scale = -5:1:5; 
        end
        v_title = 'Meridional velocity (m/s)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '04']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '15']);
        end
    elseif (strcmp(itype,'geopot')) % Plot geopotential data
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = 2400:100:3400;
        elseif (strcmp(test,'DCMIP2008c12'))
            c_scale = 750:100:1500;
        end
        v_title = 'Geopotential (m)';
        s_ll = s_ll+load([file_base itime '05']);
    elseif (strcmp(itype,'vort')) % Plot relative vorticity data
        s_ll = s_ll+load([file_base itime '06']);
        %c_scale = -3e-5:1e-5:6e-5;
        c_scale = -1e-5:1e-6:2e-5;
        v_title = 'Relative vorticity';
    elseif (strcmp(itype,'surf_press')) % Plot surface pressure data
        s_ll = s_ll+load([file_base itime '07']);
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = 800:20:1030;
        elseif (strcmp(test,'DCMIP2008c12'))
            c_scale = 930:10:1030;
        end
        v_title = 'Surface pressure';
    elseif (strcmp(itype,'eddy_mom')) % Plot zonal eddy momentum flux
        c_scale = -100:20:100; % Held-Suarez
        v_title = 'Eddy momentum flux (m^2/s^2)';
        s_zo = s_zo+load([file_base itime '16']);
    elseif (strcmp(itype,'eddy_ke')) % Plot zonal eddy kinetic energy
        if (strcmp(test,'DCMIP2008c5'))
            c_scale = 0:5:100;
        elseif (strcmp(test,'Held_Suares'))
            c_scale = 0:40:480;
        end
        v_title = 'Eddy kinetic energy (m^2/s^2)';
        s_zo = s_zo+load([file_base itime '17']);
    elseif (strcmp(itype,'eddy_heat_flux')) % Plot zonal eddy heat flux
        c_scale = -20:5:20; % Held-Suarez
        v_title = 'Eddy kinetic energy (K m/s)';
        s_zo = s_zo+load([file_base itime '18']);
    end
    % Erase extracted files
    file_erase = ['\rm ' file_base '*'];
    system(file_erase);
end

% Plot data
if (lon_lat)
    s_ll = s_ll/(t2-t1+1); % Average
    fprintf('Minimum value of variable %s = %8.4e\n',itype, min(min(s_ll)));
    fprintf('Maximum value of variable %s = %8.4e\n',itype, max(max(s_ll)));
    plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift)
    axis(ax)
end

if (zonal_avg)
    s_zo = s_zo/N; % Sample mean
    fprintf('Minimum value of variable %s = %8.4e\n',itype, min(min(s_zo)));
    fprintf('Maximum value of variable %s = %8.4e\n',itype, max(max(s_zo)));
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, smooth, 0)
end



