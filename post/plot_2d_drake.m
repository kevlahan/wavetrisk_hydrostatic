test_case = 'drake'; 
run_id    = '2layer_fill'; 

% Field to plot
set_type = 9;
if set_type == 1
    itype = 'barotropic zonal velocity';
elseif set_type == 2
    itype = 'barotropic meridional velocity'
elseif set_type == 3
    itype = 'barotropic vorticity'
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
elseif set_type == 10
    itype = 'free surface'
elseif set_type == 11
    itype = 'internal free surface'
elseif set_type == 12
    itype = 'land'
end

smooth = false;  % smooth data over two points in each direction
shift  = false;   % shift left boundary to zero longitude
lines  = false;  % plot lines

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
    c_scale = linspace(-0.6, 0.34, 100);
    v_title = 'Free surface';
elseif (strcmp(itype,'internal free surface'))
    s_ll = load([file_base '.11']);
    c_scale = linspace(-8.5, 27, 100);
    v_title = 'Internal  freesurface';
elseif (strcmp(itype,'land'))
    s_ll = load([file_base '.12']);
    c_scale = linspace(0, 1, 100);
    v_title = 'Land';
end

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
%% Plot latitude - longitude projection
figure
c_scale = linspace(-2e-5, 2e-5, 100);
plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines);
axis(ax)
system(['\rm ' run_id '.4*']);
%% Periodize in latitude
s_wrap = [s_ll(1:end-1,1:end-1); flipud(s_ll(1:end-1,1:end-1))]; 
lat_wrap = [lat(1,1:256) 180+lat(1:256)]; 
lon_wrap = lon(1:end-1);
Ly = Lx;
%% Energy spectrum: assume square array
%figure
%plot_lon_lat_data(s_wrap, lon_wrap, lat_wrap, c_scale, v_title, smooth, shift, lines);
% Energy spectrum calculation and plot

nx = size(s_wrap,1); 

n = [0 1:nx/2-1 -nx/2 -fliplr(1:nx/2-1)];
[k1, k2] = meshgrid(n, n);

kmax = nx/2;

% fft of data
fk = fft2(s_wrap)/nx;

% Energy spectrum integrated over shells
Ek = zeros(kmax+1,1);
dk = 2*pi/Lx;
for ix = 1:nx
    for iy = 1:nx
        k = round(sqrt(k1(ix,iy)^2 + k2(ix,iy)^2));
        if k <= kmax
%            Ek(k+1) = Ek(k+1) + conj(fk(ix,iy)).*fk(ix,iy);
             Ek(k+1) = Ek(k+1) + conj(fk(ix,iy)).*fk(ix,iy) ./ (dk*k)^2;
        end
    end
end

k=(0:kmax)'*dk;
loglog(k(2:end),Ek(2:end),'-r','LineWidth',1.2);hold on;grid on;
xlabel('k');ylabel('E(k)');
%%
loglog(k(4:end),k(4:end).^(-4)/4e9,'b--','LineWidth',1.2);
loglog(k(4:end),k(4:end).^(-2)/1e8,'g--','LineWidth',1.2);
loglog(k(4:end),k(4:end).^(-2)/1e7,'r--','LineWidth',1.2);
%%
legend('barotropic','layer 1 baroclinic','layer 2 baroclinic', 'k^{-4}','k^{-2}','k^{-2}');
title('Energy spectrum');set(gca,'FontSize',16);



%% 2d cwt energy spectrum
scales = 2.^(linspace(1,10,32));
cwt2d = cwtft2(s_wrap,'wavelet',{'isomorl',{6,1}},'scales',scales,'angles',1,'norm','L2');

% Average over space
Ew = zeros(numel(scales),1);
for s = 1:numel(scales)
    Ew(s) = sum(abs(cwt2d.cfs(:,:,1,s,1)).^2,[1,2]);
end

% Rescale to energy spectrum
k = 2*pi./(scales*dx)';
%Ew = Ew .* k*pi;
Ew = Ew ./k.^2 .* k*pi;
loglog(k, Ew,'g');
%% Plot cwt at a particular scale
j_scl = 3;
scl = scales(j_scl)
figure;imagesc(abs(cwt2d.cfs(:,:,1,j_scl,1)).^2);

