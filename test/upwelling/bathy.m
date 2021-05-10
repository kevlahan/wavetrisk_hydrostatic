% Bathymetry profile for upwelling case on the sphere
clear
slope = 1.8e-4;
shift = 10;
b_max = 125;

Ly = 80e3;
radius = 240e3;
lat_c = 45.2; % centre of zonal channel
lat_width = rad2deg(Ly/radius);
lat = linspace(-90,90,256);
y = (lat-lat_c+lat_width/2)/180 * pi*radius;
dy = y(2)-y(1);

A = b_max/(1-tanh(-slope*Ly/shift));

for i = 1:numel(y)
    h(i) = A * (1 - tanh(slope*(f(y(i),Ly) - Ly/shift)));
end
figure
plot(lat, h,'b-','linewidth',2); axis([45-lat_width/2 45+lat_width/2 0 150]);
xlabel('Latitude');ylabel('bathymetry (m)');grid on;axis('square')
%%
plot(lat, 125/2 * (1-cos(2*pi*(lat-lat_c)/lat_width)),'o');axis([45-lat_width/2 45+lat_width/2 0 150]);
%% ROMS vertical levels profile
z_edge = [0 1.6 6.0 8.7 10.6 12.1 13.4 14.3 15.2 16.0 16.7 17.3 17.9 18.4 19.0 19.6 19.8]-19.8;
b_vert = -z_edge'/19.8;

k = ((1:numel(b_vert))-1)/(numel(b_vert)-1);
pb = polyfit(k,b_vert,4);
z = linspace(k(1),k(end),300);
plot(k,b_vert,'o');hold on; plot(z,polyval(pb,z),'r');

%% ROMS temperature profile
z_node = 0.5*(b_vert(1:end-1) + b_vert(2:end))*(-150);

temp = 14:0.5:21.5;

rho0  = 1027;
T0    = 14;
beta  = 0.28;
rho =rho0 - beta*(temp'-T0);
plot(temp,z_node,'o');
ylabel('z');xlabel('temperature');hold on;
p = polyfit(z_node,temp,3);
z = linspace(min(z_node),max(z_node),300);
temp_fit = polyval(p,z);

rho =rho0 - beta*(temp_fit-T0);
plot(temp_fit,z,'r','linewidth',2); grid on;
figure;plot(z, rho,'b','linewidth',2); grid on;
%% CROCO profile
z_edge = [-150 -104 -74 -54 -40 -31 -24 -20 -16 -13 -11 -8.6 -6.7 -5 -3.25 -1.6 0]; % croco
z_edge = [-150 -80  -60 -55 -52 -45 -40 -35 -30 -25 -20 -15  -11 -7.5 -4.50 -2.0 0];

b_vert = -z_edge/150
k = ((1:numel(b_vert))-1)/(numel(b_vert)-1);
pb = polyfit(k,b_vert,7);
H=-150;
z = linspace(k(1),k(end),300);
b_fit = polyval(pb,z); b_fit(1)=1;b_fit(end) = 0;
plot(k,b_vert,'o');hold on; plot(z,b_fit,'r');

%% Temperature profile
T0    = 14;
H0    = -150;
hz    = 6.5;

rho0  = 1027;
beta  = 0.28;
H     = -150;

hz    = 6.5 * abs(H/150);
z0    = -35 * abs(H/150);
z1    = -75 * abs(H/150);
strat = abs(H);

z=linspace(H,0,500);
T = @(z) T0 + 4*tanh((z-z0)/hz) + (z-z1)/strat;
figure(1);plot(z,T(z),'linewidth',2);grid on;xlabel('z');ylabel('T');
set(gca,'fontsize',18);

rho = @(z) rho0 - beta*(T(z)-T0);
%rho = @(z) rho0 * (1 - beta*(T(z)-T0));
figure(2);plot(rho(z),z,'linewidth',2);xlabel('\rho(z)');ylabel('z');grid on;
set(gca,'fontsize',18);

figure(3);
temp=linspace(14,20,200);
%plot(temp,rho0+beta*temp,'linewidth',2);grid on;xlabel('T');ylabel('\rho');
plot(temp,rho0-0.28*(temp-T0),'linewidth',2);grid on;xlabel('T');ylabel('\rho');
set(gca,'fontsize',18);
%% Stretched grid

zlevels = 16;
eta = 10;
H   = -150;

% Chebyshev
b_vert(1) = 1;
for k = 2:zlevels
%b_vert(k) = (1 + cos((2*(k-1)-1)/(2*(zlevels-1)) * pi)) / 2; % chebyshev
%b_vert(k) = 1 - (k-1)/zlevels; % uniform
b_vert(k) = (1 + sin((2*(k-1)-1)/(1*(zlevels-1)) * pi)) / 2
end

b_vert(1) = 1;
b_vert(2) = -103.935/H;
b_vert(3) =  -73.66/H;
b_vert(4) =  -53.57/H;
b_vert(5) =  -40.06/H;
b_vert(6) =  -30.80/H;
b_vert(7) =  -24.28/H;
b_vert(8) =  -19.54/H;
b_vert(9) =  -15.94/H;
b_vert(10) = -13.07/H;
b_vert(11) = -10.68/H;
b_vert(12) =  -8.60/H;
b_vert(13) =  -6.71/H;
b_vert(14) =  -4.95/H;
b_vert(15) =  -3.26/H;
b_vert(16) =  -1.62/H;
b_vert(17) = 0;
a_vert = 1-b_vert;

a_vert_mass = a_vert(2:zlevels+1) - a_vert(1:zlevels);
b_vert_mass = b_vert(2:zlevels+1) - b_vert(1:zlevels);

disp(' ' )
z = a_vert*eta + b_vert*H;
disp([a_vert' b_vert' z'])
disp(' ' )
dz = a_vert_mass*eta + b_vert_mass*H;
disp([a_vert_mass' b_vert_mass' dz'])

z = a_vert*eta + b_vert*H;
plot(z,0*numel(z),'bo');grid on;
axis([H eta -1 1])

function [x] = f(y, Ly)
if y <= Ly/2
    x = y;
else
    x = Ly - y;
end
end





