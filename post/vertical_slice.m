function vertical_slice (test_case, run_id, machine, directory, time, figure_type, field)
%
% Plots vertical slice data
%
% Usage: vertical_slice('upwelling','implicit','niagara.computecanada.ca','~/proj/upwelling/J6J8_nodiff', 4, 1,'temperature')
%
% Input variables:
% test_case =   test case of data to plot (e.g. 'upwelling')
% run_id      = run id of data to plot (e.g. 'upwelling_test')
% machine     = remote machine where data is stored (e.g. 'if.mcmaster.ca',
%               'niagara.computecanada.ca')
% directory   = directory where data is stored on remote machine
% time        = time index of data to plot
% figure_type = 1 (single field), 2 (all fields)
% field       = single field to plot ('density', 'meridional', 'zonal',
%               'temperature', 'vertical')

if strcmp (test_case, 'upwelling')
    ref_density = 1027;
    radius    = 240; % radius of planet in km
    lat_c     = 45;  % centre of channel in degrees
    lat_w     = 80;  % width of zonal channel in km
    max_depth = -150;
elseif strcmp (test_case, 'jet')
    ref_density = 1027.75;
    radius    = 1000;
    lat_w     = radius/4;
    
    lat_c     = 30;
    max_depth = -4000;
elseif strcmp (test_case, 'seamount')
    ref_density = 1000;
    radius    = 6371.229/41.75;
    lat_w     = radius/4;
    
    lat_c     = 43.29;
    lon_c     = 0;
    max_depth = -5000;
end
lat_w_deg = lat_w/radius * 180/pi;

% Transfer data
file_name   = [run_id '.5' '.' sprintf( '%04d', time)];
remote_file = [directory '/' file_name  '.tgz'];
local_dir   = ['~/hydro/' test_case];

scp_cmd     = ['scp ' machine ':' remote_file ' ' local_dir '/' file_name '.tgz'];
unix (sprintf(scp_cmd));

file_tar = ['tar ' 'xf ' local_dir '/' file_name '.tgz' ' -C ' local_dir '/' ];
disp(['Uncompressing file ' file_name '.tgz']);
system(file_tar);
disp(scp_cmd)

% Load coordinates
lon = fread(fopen([local_dir '/' run_id '.5.50']),'double'); Nlon = numel(lon);
lat = fread(fopen([local_dir '/' run_id '.5.51']),'double'); Nlat = numel(lat);

xlat = fread(fopen([local_dir '/' run_id '.5.52']),[Nlat,2],'double'); %xlat = [xlat ; -xlat(1,:)];
xlon = fread(fopen([local_dir '/' run_id '.5.53']),[Nlon,2],'double');

zlat = fread(fopen([local_dir '/' run_id '.5.54']),'double'); zlat = reshape(zlat,Nlat,[],2);
zlon = fread(fopen([local_dir '/' run_id '.5.55']),'double'); zlon = reshape(zlon,Nlon,[],2);

lat_slice = fread(fopen([local_dir '/' run_id '.5.56']),'double'); lat_slice = reshape(lat_slice,Nlat,[],5);
lon_slice = fread(fopen([local_dir '/' run_id '.5.57']),'double'); lon_slice = reshape(lon_slice,Nlon,[],5);

% Plot results
if figure_type == 1
    ax1=figure;plot_field (xlat, zlat, lat_slice, field, 'raw', ax1)
elseif figure_type == 2
    str = "Upwelling results at day "+compose("%1.0f",time/2)+" for J7 (dx = 2.26 km, dt = 1840 s)";
    %sgtitle(str);
    type = 'raw'; %'raw' or 'interp'
    
    figure
    ax1 = subplot(2,2,1); plot_field (xlat, zlat, lat_slice, 'temperature', type, ax1)
    ax2 = subplot(2,2,2); plot_field (xlat, zlat, lat_slice, 'zonal',       type, ax2)
    ax3 = subplot(2,2,3); plot_field (xlat, zlat, lat_slice, 'meridional',  type, ax3)
    ax4 = subplot(2,2,4); plot_field (xlat, zlat, lat_slice, 'vertical',    type, ax4)
elseif figure_type == 3
    % Plot vertical grid
    skip = 1;
    plot(lat(1:skip:end),zlat(1:skip:end,:,1)','k-','linewidth',1.2); hold on;
    plot(lat(1:skip:end),zlat(1:skip:end,:,2)','k-','linewidth',1.2);
    axis([lat_c-lat_w_deg/2 lat_c+lat_w_deg/2 max_depth 0]);
    plot([lat(1:skip:end) lat(1:skip:end)], ylim, 'k-','linewidth',1.2);
    %axis([0 90 -150 0])
    axis([min(xlat(:,2,1)) max(xlat(:,1,1)) max_depth 0])
    xlabel('latitude'); ylabel('z (m)'); hold on;
    set(gca,'fontsize',18);hold off
end
    function plot_field (xlat, zlat, lat_slice, field, type, ax)
        nz = 12; % interpolate to this number of vertical levels
        skip = 1;
        zlevels = size(zlat,2);
        Nlat    = size(zlat,1);
        DAY = 60^2 * 24;
        if strcmp(field,'density') % used for seamount
            m = 3;
            if strcmp (test_case, 'upwelling')
                c1 = linspace(-3, 0, 20);
                c2 = linspace(-3, 0, 20);
                trans = @(rho) rho - 1027;
            elseif strcmp (test_case, 'jet')
                c1 = linspace(24, 28, 200);
                c2 = 24:0.5:28;
                trans = @(rho) rho - 1000;
            elseif strcmp (test_case, 'seamount')
                c1 = linspace(-3, 0, 200);
                c2 = linspace(-3, 0, 20);
                trans = @(rho) rho - 1000;
            end
        elseif strcmp(field,'temperature')
            m = 3;
            c1 = linspace(9, 19, 200);
            c2 = [10 11 12 13 14 15 16 17 18 19 20];
            trans = @(rho) 14 + (1027-rho)/0.28;
        elseif strcmp(field,'zonal')
            m = 1;
            if strcmp (test_case, 'upwelling')
                c1 = linspace(-20, 20, 200);
                c2 = [-18 -16 -14 -12 -10 -8 -6 -4 -2 0];
                trans = @(u) -100 * u;
            elseif strcmp (test_case, 'jet')
                c1 = linspace(-70, 70, 200);
                c2 = [5 15 25 35 45 55];
                trans = @(u) 100 * u;
            end
        elseif strcmp(field,'meridional')
            m = 2;
            c1 = linspace(-6, 6, 200);
            c2 = [-6 -5 -4 -3 -2 -1  1 2 3 4];
            trans = @(u) 100 * u;
        elseif strcmp(field,'vertical')
            m = 5;
            c1 = linspace(-15, 15, 200);
            c2 = [-14 -12 -10 -8 -6 -4 -2  2 4 6 8 10 12 14];
            trans = @(u) DAY * u;
        end
        
        dat = trans(lat_slice(:,:,m));
                
        %c = linspace(min(dat(:)), max(dat(:)), 100);
        
        if strcmp(type,'interp')
            % Croco grid
            p = [-5.7831  18.9754 -24.6521  16.1698 -5.7092 0.9972];
            b_vert(1) = 1; b_vert(nz+1) = 0;
            for k = 2:nz
                b_vert(k) = polyval(p,(k-1)/nz);
            end
            %b_vert_mass = 0.5 * (b_vert(2:nz+1) + b_vert(1:nz));
            
            x_node = repelem(0.5*(xlat(:,1)+xlat(:,2)),1,zlevels+2);
            z_node = [zlat(:,1,1) 0.5*(zlat(:,:,1)+zlat(:,:,2)) zlat(:,zlevels,2)];
            x_unif = repelem(0.5*(xlat(:,1)+xlat(:,2)),1,nz+1); x_unif = x_unif(1:skip:end,:);
            
            for i = 1:size(x_unif,1)
                j = 1 + skip*(i-1);
                %         z_unif(i,:) = linspace(min(z_node(j,:)),max(z_node(j,:)), nz);
                z_s = min(z_node(j,:));
                for k = 1:nz+1
                    z_unif(i,k) = b_vert(k) * z_s;
                end
            end
            %dat = smooth2a([dat(:,1) dat dat(:,zlevels)],2,2);
            dat = [dat(:,1) dat dat(:,zlevels)];
            data = smooth2a(griddata(x_node, z_node, dat, x_unif, z_unif),2,2);
            contourf(x_unif, z_unif, data, c1, 'LineColor', 'none'); hold on
        elseif strcmp(type,'raw')
            for k = 1:zlevels
                z11 = zlat(1,k,1);
                z12 = 0.5 * (zlat(1,k,1) + zlat(2,k,1));
                z21 = 0.5 * (zlat(1,k,2) + zlat(2,k,2));
                z22 = zlat(1,k,2);
                x = [xlat(1,1) xlat(1,2) xlat(1,2) xlat(1,1)];
                y = [z11 z12 z21 z22];
                patch(x,y,dat(1,k),'LineStyle','none');
                for i = 2:Nlat-1
                    z11 = 0.5 * (zlat(i-1,k,1) + zlat(i,  k,1));
                    z12 = 0.5 * (zlat(i,  k,1) + zlat(i+1,k,1));
                    
                    z21 = 0.5 * (zlat(i,  k,2) + zlat(i+1,k,2));
                    z22 = 0.5 * (zlat(i-1,k,2) + zlat(i,  k,2));
                    
                    x = [xlat(i,1) xlat(i,2) xlat(i,2) xlat(i,1)];
                    y = [z11 z12 z21 z22];
                    patch(x,y,dat(i,k),'LineStyle','none');
                end
                z11 = 0.5 * (zlat(Nlat-1,k,1) + zlat(Nlat,k,1));
                z12 = zlat(Nlat,k,1);
                z21 = zlat(Nlat,k,2);
                z22 = 0.5 * (zlat(Nlat-1,k,2) + zlat(Nlat,k,2));
                x = [xlat(Nlat,1) xlat(Nlat,2) xlat(Nlat,2) xlat(Nlat,1)];
                y = [z11 z12 z21 z22];
                patch(x,y,dat(Nlat,k));
            end
            hold on
        end
        
        if strcmp(field,'temperature')
            colormap(ax, jet)
        elseif strcmp(field,'density')
            colormap(ax, autumn)
        else
            colormap(ax, hsv)
        end
        
        % Contour lines from raw data
        z = 0.5*(zlat(:,:,1) + zlat(:,:,2));
        x = repelem(0.5*(xlat(:,1)+xlat(:,2)),1,size(z,2));
        dat = trans(lat_slice(:,:,m));
        %contour(x, z, dat, c2, 'k-',  'Linewidth', 1.5);
        %contour(x, z, dat, min(c2) : max(c2(c2<0)),   'k--', 'Linewidth', 1.5);
        %contour(x, z, dat, min(c2(c2>=0))  : max(c2), 'k-',  'Linewidth', 1.5);
        
        max(xlat(:,2,1))
        
        xlabel('latitude'); ylabel('z (m)');
        
        caxis([min(c1) max(c1)]);
        axis([min(xlat(:,2,1)) max(xlat(:,1,1)) max_depth 0])
        %axis([lat_c-lat_w_deg/2 lat_c+lat_w_deg/2 max_depth 0]);
        set(gca,'fontsize',18);
        
        hcb=colorbar;
        if strcmp(field,'density')
            hcb.Label.String = "Density (kg/m^3)";
        elseif strcmp(field,'temperature')
            hcb.Label.String = "Temperature (celsius)";
        elseif strcmp(field,'zonal')
            hcb.Label.String = "Zonal velocity (cm/s)";
        elseif strcmp(field,'meridional')
            hcb.Label.String = "Meridional velocity (cm/s)";
        elseif strcmp(field,'vertical')
            hcb.Label.String = "Vertical velocity (m/day)";
        end
        
        fprintf('Minimum value of %s = %8.4e\n', field, min(min(dat)));
        fprintf('Maximum value of %s = %8.4e\n', field, max(max(dat)));
    end
end
