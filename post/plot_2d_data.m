% Plot 2d data from export_2d
clear; %close all;

% Test case
test_case = 'DCMIP2012c4';
%test_case = 'DCMIP2008c5';
%test_case = 'Held_Suarez';

machine   = 'if';
%machine   = 'mac';
t1        = 9; % Start time
t2        = t1; % End time
% Options: 'temp' 'zonal' 'merid' 'geopot' 'vort' 'surf_press' 'temp_var' 'eddy_mom' 'eddy_ke' 'eddy_heat_flux'
itype     = 'vort';

lon_lat   = true; % Plot longitude - latitude data
zonal_avg = false; % Plot zonally averaged data
shift     = true; % shift left boundary to zero longitude
smooth    = false; % smooth data over two points in each direction
lines     = false; % remove lines

if (strcmp(itype,'temp_var')||strcmp(itype,'eddy_mom')||strcmp(itype,'eddy_ke')||strcmp(itype,'eddy_heat_flux'))
    lon_lat = 0;
    zonal_avg = 1;
end

% Axis limits
if (strcmp(test_case,'DCMIP2008c5'))
    ax = [0 360 -90 90];
elseif (strcmp(test_case,'DCMIP2012c4'))
    if (strcmp(itype,'temp')||strcmp(itype,'surf_press'))
        ax = [45 360 0 90];
    else
        ax = [120 270 25 75];
        %ax = [90 200 25 75];
        %ax = [0 360 25 75];
    end
end

N = t2-t1+1; % number of samples

file_base = [test_case '.3'];
if (strcmp(machine,'if'))
    pathid = ['/net/if/1/home/kevlahan/data/jobs/' test_case '/'];
elseif (strcmp(machine,'mac'))
    pathid = ['/Users/kevlahan/hydro/' test_case '/'];
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
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = 270:3:303;
        elseif (strcmp(test_case,'DCMIP2012c4'))
            c_scale = 220:10:310;
        elseif (strcmp(test_case,'Held_Suarez'))
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
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = -15:5:50; 
        elseif (strcmp(test_case,'DCMIP2012c4'))
            c_scale = -60:5:60;   
        elseif (strcmp(test_case,'Held_Suarez'))
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
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = -35:5:20;
        elseif (strcmp(test_case,'DCMIP2012c4'))
            c_scale = -35:5:35; 
        end
        v_title = 'Meridional velocity (m/s)';
        if (lon_lat)
            s_ll = s_ll+load([file_base itime '04']);
        end
        if (zonal_avg)
            s_zo = s_zo+load([file_base itime '15']);
        end
    elseif (strcmp(itype,'geopot')) % Plot geopotential data
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = 2400:100:3400;
        elseif (strcmp(test_case,'DCMIP2012c4'))
            c_scale = 750:100:1500;
        end
        v_title = 'Geopotential (m)';
        s_ll = s_ll+load([file_base itime '05']);
    elseif (strcmp(itype,'vort')) % Plot relative vorticity data
        s_ll = s_ll+load([file_base itime '06']);
         if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = -3e-5:1e-5:6e-5;
        elseif (strcmp(test_case,'DCMIP2012c4'))
            %if (t<=48)
            %    c_scale = -3e-5:1e-5:6e-5;
            %    ax = [90 200 25 75];
            %else
                c_scale = -10e-5:5e-5:40e-5;
                ax = [120 270 25 75];
            %end
        end
        v_title = 'Relative vorticity';
    elseif (strcmp(itype,'surf_press')) % Plot surface pressure data
        s_ll = s_ll+load([file_base itime '07']);
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = 800:20:1030;
        elseif (strcmp(test_case,'DCMIP2012c4'))
            c_scale = 930:10:1030;
        end
        v_title = 'Surface pressure';
    elseif (strcmp(itype,'eddy_mom')) % Plot zonal eddy momentum flux
        c_scale = -100:20:100; % Held-Suarez
        v_title = 'Eddy momentum flux (m^2/s^2)';
        s_zo = s_zo+load([file_base itime '16']);
    elseif (strcmp(itype,'eddy_ke')) % Plot zonal eddy kinetic energy
        if (strcmp(test_case,'DCMIP2008c5'))
            c_scale = 0:5:100;
        elseif (strcmp(test_case,'Held_Suares'))
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
    fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_ll)));
    fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_ll)));
    figure;plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines)
    axis(ax)
end

if (zonal_avg)
    s_zo = s_zo/N; % Sample mean
    fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_zo)));
    fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_zo)));
    plot_zonal_avg_data(s_zo, lat, P_z, c_scale, v_title, smooth, 0)
end



