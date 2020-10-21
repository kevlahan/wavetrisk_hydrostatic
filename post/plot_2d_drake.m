%% Load file
N         = 2048;  % number of longitude points in projection
test_case = 'drake';
run_id    = '1layer_J6'; zlevels   = 1;
% run_id  = '2layer_J6'; zlevels   = 2;
itype   = 'barotropic vorticity' % field to analyze
machine = 'niagara.computecanada.ca'

% itype values:
% 'barotropic zonal velocity' .    (zonal velocity for single layer case)
% 'barotropic meridional velocity' (meridional velocity for single layer case)
% 'barotropic vorticity'           (vorticity for single layer case)
% 'layer 1 baroclinic zonal velocity'
% 'layer 1 baroclinic meridional velocity'
% 'layer 1 baroclinic vorticity'
% 'layer 2 baroclinic zonal velocity'
% 'layer 2 baroclinic meridional velocity'
% 'layer 2 baroclinic vorticity'
% 'free surface'
% 'internal free surface'
% 'land'

file_base   = [run_id '.4'];
remote_file = ['~/hydro/' test_case '/' file_base '.tgz'];
local_file  = ['~/hydro/' test_case '/' file_base '.tgz'];
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];

unix (sprintf(scp_cmd));

file_tar = ['tar ' 'xf ' local_file];
disp(['Uncompressing file ' local_file]);
system(file_tar);

% Load coordinates
lon = fread(fopen([file_base '.20']),'double');
lat = fread(fopen([file_base '.21']),'double');

% Load field to analyze
smooth = false;  % smooth data over two points in each direction
shift  = false;   % shift left boundary to zero longitude
lines  = false;  % plot lines

if (strcmp(itype,'barotropic zonal velocity'))
    fid = fopen([file_base '.01']);
    c_scale = linspace(-0.8, 1.5, 100);
    v_title = 'Barotropic zonal velocity';
elseif (strcmp(itype,'barotropic meridional velocity'))
    fid = fopen([file_base '.02']);
    c_scale = linspace(-1.1, 1.1, 100);
    v_title = 'Barotropic meridional velocity';
elseif (strcmp(itype,'barotropic vorticity'))
    fid = fopen([file_base '.03']);
    c_scale = linspace(-3e-5, 3e-5, 100);
    v_title = 'Barotropic vorticity';
elseif (strcmp(itype,'layer 1 baroclinic zonal velocity'))
    fid = fopen([file_base '.04']);
    c_scale = linspace(-0.07, 0.04, 100);
    v_title = 'Baroclinic zonal velocity';
elseif (strcmp(itype,'layer 1 baroclinic meridional velocity'))
    fid = fopen([file_base '.05']);
    c_scale = linspace(-0.06, 0.05, 100);
    v_title = 'Baroclinic meridional velocity';
elseif (strcmp(itype,'layer 1 baroclinic vorticity'))
    fid = fopen([file_base '.06']);
    c_scale = linspace(-5e-6, 5e-6, 100);
    v_title = 'Layer 1 baroclinic vorticity';
elseif (strcmp(itype,'layer 2 baroclinic zonal velocity'))
    fid = fopen([file_base '.07']);
    c_scale = linspace(-0.1, 0.2, 100);
    v_title = 'Baroclinic zonal velocity';
elseif (strcmp(itype,'layer 2 baroclinic meridional velocity'))
    fid = fopen([file_base '.08']);
    c_scale = linspace(-0.13, 0.16, 100);
    v_title = 'Baroclinic meridional velocity';
elseif (strcmp(itype,'layer 2 baroclinic vorticity'))
    fid = fopen([file_base '.09']);
    c_scale = linspace(-2e-5, 2e-5, 100);
    v_title = 'Layer 2 baroclinic vorticity';
elseif (strcmp(itype,'free surface'))
    if zlevels == 2
        fid = fopen([file_base '.10']);
    elseif zlevels == 1
        fid = fopen([file_base '.04']);
    end
    c_scale = linspace(-0.6, 0.34, 100);
    v_title = 'Free surface';
elseif (strcmp(itype,'internal free surface'))
    fid = fopen([file_base '.11']);
    c_scale = linspace(-8.5, 27, 100);
    v_title = 'Internal  freesurface';
elseif (strcmp(itype,'land'))
    if zlevels == 2
        fid = fopen([file_base '.12']);
    elseif zlevels == 1
        fid = fopen([file_base '.05']);
    end
    c_scale = linspace(0, 1, 100);
    v_title = 'Land';
end
s_ll = reshape (fread(fid,'double'), N+1, N/2+1)';

N         = size(s_ll,2); % resolution of projection
radius    = 6371.229e3/6; % radius of Earth in metres

Nx        = N;     % number of points in longitude
Ny        = N/2;   % number of points in latitude
Lx        = 2*pi*radius; 
Ly        = pi*radius; 
dx        = Lx/Nx;
dy        = Ly/Ny;

lat_min   = min(lat);
lat_max   = max(lat);
lon_min   = min(lon);
lon_max   = max(lon);
ax = [lon_min lon_max lat_min lat_max];

fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_ll)));
fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_ll)));

% Erase files
file_erase = ['\rm ' file_base '*'];
system(file_erase);
%% Plot latitude - longitude projection
figure
c_scale = linspace(-3.5e-5, 3.5e-5, 100);
plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines);
axis(ax)
%%
%print -dpng ~/hydro/drake/barotropic_proj.png
%print -dpng ~/hydro/drake/barotropic1.png
print -dpng ~/hydro/drake/barotropic2.png
%% Energy spectrum: assume square array

% Periodize in latitude
s_wrap = [s_ll(1:end-1,1:end-1); flipud(s_ll(1:end-1,1:end-1))]; 
lat_wrap = [lat(1:size(s_ll,1)-1) 180+lat(1:size(s_ll,1)-1)]; 
lon_wrap = lon(1:end-1);
Ly = Lx;

% Plot wrapped data
%figure
%plot_lon_lat_data(s_wrap, lon_wrap, lat_wrap, c_scale, v_title, smooth, shift, lines);

nx = size(s_wrap,1); 
n = [0 1:nx/2-1 -nx/2 -fliplr(1:nx/2-1)];
[k1, k2] = meshgrid(n, n);

kmax = nx/2;
dk   = 1; % assume L = 2*pi
k    = (0:kmax)'*dk;
Ek   = zeros(kmax+1,1); % initialize energy spectrum

% Energy spectrum integrated over shells
fk = fft2(s_wrap)/nx; % fft of data
for ix = 1:nx
    for iy = 1:nx
        K = round(sqrt(k1(ix,iy)^2 + k2(ix,iy)^2));
        if K <= kmax
            %Ek(K+1) = Ek(K+1) + 0.5 * conj(fk(ix,iy)).*fk(ix,iy);
            Ek(K+1) = Ek(K+1) + 0.5 * conj(fk(ix,iy)).*fk(ix,iy) * (dx./(dk*K))^2;
        end
    end
end

% Check energy conservation
%sum(s_wrap.^2,[1,2]) - sum(Ek)

if strcmp (itype, 'barotropic vorticity')
    col = 'b';
elseif strcmp (itype, 'layer 1 baroclinic vorticity')
    col = 'm';
elseif strcmp (itype, 'layer 2 baroclinic vorticity')
    col = 'r';
end

% Energy spectrum
loglog(k(2:end),Ek(2:end),col,'LineWidth',1.2,'DisplayName',itype); hold on;

% Power law scalings
if (strcmp(itype,'layer 2 baroclinic vorticity'))
    k1 = 4; k2 = round(numel(k)/8);
    loglog(k(k1:k2),k(4:k2).^(-5/3)/1.5e1,'b--','LineWidth',1.2,'DisplayName','k^{-5/3}');
    k1 = 60; k2 = round(numel(k)*0.9);
    loglog(k(k1:k2),k(60:k2).^(-4)*4e3,'g--','LineWidth',1.2,'DisplayName','k^{-4}');
end

if (strcmp(itype,'layer 1 baroclinic vorticity'))
    k1 = 4; k2 = round(numel(k)/8);
    loglog(k(k1:k2),k(4:k2).^(-5/3)/1.5e2,'b--','LineWidth',1.2,'DisplayName','k^{-5/3}');
    k1 = 60; k2 = round(numel(k)*0.9);
    loglog(k(k1:k2),k(60:k2).^(-4)*4e2,'g--','LineWidth',1.2,'DisplayName','k^{-4}');
end

if (strcmp(itype,'barotropic vorticity'))
    if zlevels == 2
        k1 = 4; k2 = round(numel(k)/4);
        loglog(k(k1:k2),k(k1:k2).^(-3.3)*2e2,'r--','LineWidth',1.2,'DisplayName','k^{-3.3}');
        k1 = 100; k2 = round(numel(k));
        loglog(k(k1:k2),k(k1:k2).^(-5)*1e6,'g--','LineWidth',1.2,'DisplayName','k^{-5}');
    elseif zlevels == 1
        k1 = 4; k2 = round(numel(k)/4);
        loglog(k(k1:k2),k(k1:k2).^(-3)*5e1,'r--','LineWidth',1.2,'DisplayName','k^{-3}');
        k1 = 100; k2 = round(numel(k));
        loglog(k(k1:k2),k(k1:k2).^(-5)*2e6,'g--','LineWidth',1.2,'DisplayName','k^{-5}');
    end
end

legend

xlabel('k');ylabel('E(k)');grid on;
title('Energy spectrum');set(gca,'FontSize',16);
%%
Ek_av = [k Ek];
save('Ek_layer2_baroclinic_av.mat','Ek_av');
%%
print -dpng ~/hydro/drake/energy_spectrum.png
%% CWT for a subset of Earth
selection = 'Laminar'
if strcmp (selection, 'Vortical')
    lon_min = -40;
    lon_max =  80;
    lat_min = -50;
    lat_max = -15;
elseif strcmp (selection, 'Filament')
    lon_min =  90;
    lon_max = 270;
    lat_min = -90;
    lat_max =  90;
elseif strcmp (selection, 'Laminar')
    lon_min = -150;
    lon_max = -100;
    lat_min = -90;
    lat_max =  90;
elseif strcmp (selection, 'Entire')
    lon_min = -180;
    lon_max =  180;
    lat_min = -90;
    lat_max =  90+180;
end

A = find(abs(lat_wrap-lat_min)<180/(nx)); coord_lat(1) = A(1)
A = find(abs(lat_wrap-lat_max)<180/(nx)); coord_lat(2) = min([A numel(lat_wrap)]);

A = find(abs(lon_wrap-lon_min)<360/nx); coord_lon(1) = A(1)
A = find(abs(lon_wrap-lon_max)<360/nx); coord_lon(2) = min([A numel(lon_wrap)]);

lat_sel = lat_wrap(coord_lat(1):coord_lat(2));
lon_sel = lon_wrap(coord_lon(1):coord_lon(2));
s_sel = s_wrap(coord_lat(1):coord_lat(2), coord_lon(1):coord_lon(2));

% Show selection
figure; 
c_scale = linspace(min(s_sel(:)), max(s_sel(:)), 100); 
colormap(jet(numel(c_scale)-1));caxis([min(c_scale) max(c_scale)]);c=colorbar;
c.Label.String=v_title;c.Label.FontSize=16;c.YTick=c_scale;
axis('equal'); 
[~,h]=contourf(lon_sel, lat_sel, s_sel, c_scale, 'LineColor','none');
xlabel('Longitude');ylabel('Latitude');
title([selection ' selection for wavelet spectrum']);set(gca,'FontSize',16);
%% 
print_save = ["~/hydro/drake/selection_" selection ".png"]
%print -dpng ~/hydro/drake/laminar_selection.png
print -dpng ~/hydro/drake/vortical_selection.png

%% 2d cwt energy spectrum
n_select = min(coord_lat(2)-coord_lat(1)+1, coord_lon(2)-coord_lon(1)+1)
scales = 2.^(linspace(1,log(n_select)/log(2),32));
cwt2d = cwtft2(s_wrap,'wavelet',{'isomorl',{6,1}},'scales',scales,'angles',1,'norm','L2');
%cwt2d = cwtft2(s_wrap,'wavelet','cauchy','scales',scales,'angles',1,'norm','L2');

% Average over space
Ew = zeros(numel(scales),1);
for s = 1:numel(scales)
    Ew(s) = sum(abs(cwt2d.cfs(coord_lat(1):coord_lat(2), coord_lon(1):coord_lon(2), 1, s, 1)).^2,[1,2]);
end
Ew = 0.5 * Ew * (dx/nx)^2;

% Rescale to energy spectrum integrated over shells
k_wav = 2*pi./(scales * 2*pi/nx)';

%Ew = Ew * pi.*k_wav;
Ew = Ew ./k_wav.^2 * pi.*k_wav;
loglog(k_wav, Ew,'g','DisplayName', [selection ' selection']); hold on;
loglog(k(2:end),Ek(2:end),col,'LineWidth',1.2,'DisplayName',itype); hold on;
legend
xlabel('k');ylabel('E(k)');grid on;
title('Energy and averaged wavelet spectra');set(gca,'FontSize',16);
%%
print -dpng ~/hydro/drake/wavelet_spectra.png
%% Plot cwt at a particular scale
j_scl = 3;
scl = scales(j_scl)
figure;colormap('jet');imagesc(abs(cwt2d.cfs(:,:,1,j_scl,1)).^2);
