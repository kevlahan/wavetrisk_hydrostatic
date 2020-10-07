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
    c_scale = linspace(-2.4e-6, 2.4e-6, 100);
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
plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines);
axis(ax)
system(['\rm ' run_id '.4*']);
%% Energy spectrum

% Periodize in latitude
s_wrap = [s_ll(1:end-1,1:end-1); flipud(s_ll(1:end-1,1:end-1))]; 
lat_wrap = [lat(1,1:256) 180+lat(1:256)]; 
lon_wrap = lon(1:end-1);
Ly = Lx;

% figure
% plot_lon_lat_data(s_wrap, lon_wrap, lat_wrap, c_scale, v_title, smooth, shift, lines);

nx = size(s_wrap,1); ny = size(s_wrap,2);

n1 = [0 1:nx/2-1 -nx/2 -fliplr(1:nx/2-1)];
n2 = [0 1:ny/2-1 -ny/2 -fliplr(1:ny/2-1)];
[n1, n2] = meshgrid(n1, n2);
 
dk1 = pi/dx * [1 0];
dk2 = pi/dy * [0 1];

kvec = zeros(nx,ny,2);
for ix=1:nx
    for iy=1:ny
        for jj=1:2
            kvec(ix,iy,jj) = n1(ix,iy)*dk1(jj) + n2(ix,iy)*dk2(jj);
        end
    end
end
nn = min(nx, ny);
kmax = min(max(max(kvec)))/dk;
dk = pi/dx;

% fft of data
fk = fft2(s_wrap) / numel(s_wrap);

% Energy spectrum integrated over shells
Ek = zeros(kmax+1,1);
for ix = 1:nx
    for iy = 1:ny
        k = round(sqrt(kvec(ix,iy,1)^2 + kvec(ix,iy,2)^2)/dk);
        if k <= kmax
            Ek(k+1) = Ek(k+1) + conj(fk(ix,iy)).*fk(ix,iy)/(k*dk)^2;
        end
    end
end
Ek = Ek/2*dx^2;
k=(1:kmax)'*dk;
loglog(k,Ek(2:end),'-r','LineWidth',1.2);hold on;grid on;
xlabel('k');ylabel('E(k)');
%%
loglog(k(4:end),k(4:end).^(-3)/1e10,'-k','LineWidth',1.2);
%loglog(k(4:end),k(4:end).^(-2)/1e8,'-k','LineWidth',1.2);
%%
legend('barotropic','layer 1 baroclinic','layer 2 baroclinic', 'k^{-3}');
title('Energy spectrum');
