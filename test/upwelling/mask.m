% Mask function for a zonal channel

lat_c      = 30; % centre of channel (latitude in degrees)
npts_penal = 2.5; % number of points to smooth over
dx_max     = 4.5; % grid size (km)
radius     = 240; % planet radius (km)
width      = 80;  % meridional extent of channel (km)

deg        = pi/180;
lat_width  = width/radius / deg;
lat_width  = 2*10;
lat        = linspace(-90,90,1000);

% Widening of channel
dlat = 0.5*npts_penal * (dx_max/radius) / deg;

% Extent of land mass from south and north polew
width_S = 90 + (lat_c - (lat_width/2 + dlat));
width_N = 90 - (lat_c + (lat_width/2 + dlat));

% Smoothing
n_smth_S = 4*radius * width_S*deg / (dx_max * npts_penal);
n_smth_N = 4*radius * width_N*deg / (dx_max * npts_penal);

% Masks for each land mass
chi_S = exp(-abs((lat+90)/width_S).^n_smth_S); 
chi_N = exp(-abs((90-lat)/width_N).^n_smth_N); 

plot(lat, chi_S,'b'); hold on;
plot(lat,lat < lat_c-lat_width/2,'r');

plot(lat, chi_N,'b'); hold on;
plot(lat,lat > lat_c+lat_width/2,'r');

xlabel('latitude');grid on;