clear all
% test_case = 'seamount';
% run_id    = 'sea_drho_3_nonadapt';

test_case = 'upwelling';
run_id    = 'upwelling_J7';
time      = 1;
radius    = 120; % radius of planet in km
lat_w     = 80;  % width of zonal channel in km
lat_w_deg = lat_w/radius * 180/pi;

machine   = 'if.mcmaster.ca';
%machine   = 'cherry';

% Transfer data
directory   = ['~/hydro/' test_case];
file_base   = [run_id '.5'];

if ~strcmp(machine,'cherry')
remote_file = [directory '/' file_base '.' sprintf( '%04d', time) '.tgz'];
local_file  = remote_file;
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
unix (sprintf(scp_cmd));
file_tar = ['tar ' 'xf ' local_file ' -C ' directory];
disp(['Uncompressing file ' local_file]);
system(file_tar);
end

% Load coordinates
lon = fread(fopen([directory '/' file_base '.50']),'double'); Nlon = numel(lon);
lat = fread(fopen([directory '/' file_base '.51']),'double'); Nlat = numel(lat); 

xlat = fread(fopen([directory '/' file_base '.52']),[Nlat,2],'double'); %xlat = [xlat ; -xlat(1,:)];
xlon = fread(fopen([directory '/' file_base '.53']),[Nlon,2],'double');

zlat = fread(fopen([directory '/' file_base '.54']),'double'); zlat = reshape(zlat,Nlat,[],2); 
zlon = fread(fopen([directory '/' file_base '.55']),'double'); zlon = reshape(zlon,Nlon,[],2); 

lat_slice = fread(fopen([directory '/' file_base '.56']),'double');lat_slice = reshape(lat_slice,Nlat,[],5);
lon_slice = fread(fopen([directory '/' file_base '.57']),'double');lon_slice = reshape(lon_slice,Nlon,[],5);

type='raw';
figure;plot_field (xlat, zlat, lat_slice, lat_w_deg, 'vertical', type)
%% Plot results
figure
sgtitle('Upwelling results at day 5 for J7 Npts = 1.5 K_m = 2e-2')
type = 'raw'; % 'raw' or 'interp'
subplot(2,2,1); plot_field (xlat, zlat, lat_slice, lat_w_deg, 'temperature', type)
subplot(2,2,2); plot_field (xlat, zlat, lat_slice, lat_w_deg, 'zonal',       type)
subplot(2,2,3); plot_field (xlat, zlat, lat_slice, lat_w_deg, 'meridional',  type)
subplot(2,2,4); plot_field (xlat, zlat, lat_slice, lat_w_deg, 'vertical',    type)

%% Plot vertical grid
plot(lat,zlat(:,:,1)',  'k-','linewidth',1.2); hold on;
plot(lat,zlat(:,end,2)','k-','linewidth',1.2); 
axis([45-lat_w_deg/2 45+lat_w_deg/2 -150 0]);
plot([lat lat], ylim, 'k-','linewidth',1.2);
%axis([0 90 -150 0])
xlabel('latitude'); ylabel('z (m)'); hold on;
set(gca,'fontsize',18);hold off

%%
function plot_field (xlat, zlat, lat_slice, lat_width, field, type)

zlevels = size(zlat,2);
Nlat    = size(zlat,1);
DAY = 60^2 * 24;
if strcmp(field,'density')
    m = 3;
    c_scale = linspace(1025.5, 1028, 200);
    trans = @(rho) rho;
elseif strcmp(field,'temperature')
    m = 3;
    c_scale = linspace(9.5, 18.5, 200); 
    trans = @(rho) 14 + (1027-rho)/0.28;
elseif strcmp(field,'zonal')
    m = 1;
    c_scale = linspace(-35, 35, 200); 
    trans = @(u) 100 * u;
elseif strcmp(field,'meridional')
    m = 2;
    c_scale = linspace(-10, 10, 200); 
    trans = @(u) 100 * u;
elseif strcmp(field,'vertical')
    m = 5;
    c_scale = linspace(-150, 150, 200); 
    trans = @(u) DAY * u;
end

dat = trans(lat_slice(:,:,m));
%c_scale = linspace(min(dat(:)), max(dat(:)), 100);

if strcmp(type,'interp')
    x_node = repelem(0.5*(xlat(:,1)+xlat(:,2)),1,zlevels+2);
    z_node = [zlat(:,1,1) 0.5*(zlat(:,:,1)+zlat(:,:,2)) zlat(:,zlevels,2)];
    
    nz = zlevels;
    skip = 2;
    x_unif = repelem(0.5*(xlat(:,1)+xlat(:,2)),1,nz); x_unif = x_unif(1:skip:end,:);
    for i = 1:size(x_unif,1)
        j = 1 + skip*(i-1);
        z_unif(i,:) = linspace(min(z_node(j,:)),max(z_node(j,:)), nz);
    end
%     dat = smooth2a([dat(:,1) dat dat(:,zlevels)],2,2);
    dat = [dat(:,1) dat dat(:,zlevels)];
    data = griddata(x_node, z_node, dat, x_unif, z_unif);
    contourf(x_unif, z_unif, data, c_scale, 'LineColor', 'none' );
elseif strcmp(type,'raw')
    for i = 1:Nlat
        for k = 1:zlevels
            x = [xlat(i,1) xlat(i,2) xlat(i,2) xlat(i,1)];
            y = [zlat(i,k,1) zlat(i,k,1) zlat(i,k,2) zlat(i,k,2)];
            patch(x,y,dat(i,k),'LineStyle','none');
        end
    end
end

xlabel('latitude'); ylabel('z (m)'); colormap(flipud(jet))
if strcmp(field,'density')
    colormap(flipud(jet(numel(c_scale)-1)));
else
    colormap(jet(numel(c_scale)-1));
end
caxis([min(c_scale) max(c_scale)]);
axis([45-lat_width/2 45+lat_width/2 -150 0]);
%axis([0 90 -155 10])
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

