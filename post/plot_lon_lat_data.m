% Plot longitude-latitude data
clear all; close all;
itime      = '000';
unif_grid  = false;
n_contours = 15;

% Extract files
file_base = 'fort.3';
file_tar = ['tar ' 'xvf ' '/Users/kevlahan/hydro/' file_base itime '.tgz'];
system(file_tar);

% Load coordinates
lon = load([file_base itime '21']);
lat = load([file_base itime '22']);
P_z = load([file_base itime '23']); % Pressure-based vertical coordinates

% Plot mass density data
s = load([file_base itime '01']); 
figure(1);colormap(jet(n_contours))
contourf(lon,lat,s,n_contours);hold on; 
c=colorbar;c.Label.String='Mass density';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight'); 
xlabel('Longitude','fontsize',16);
ylabel('Latitude','fontsize',16);

% Plot zonally averaged mass density
s = load([file_base itime '11']);
figure(2); colormap(jet(n_contours))
if (unif_grid)
    %Interpolate onto uniform pressure grid
    P_unif = linspace(min(P_z),max(P_z),numel(P_z));
    [lat_unif,P_unif] = meshgrid(lat,P_unif);
    contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), n_contours);
    c=colorbar;c.Label.String='Mass density';c.Label.FontSize=12;
else
    contourf(lat,P_z,s);
    c=colorbar;c.Label.String='Mass density';c.Label.FontSize=12;
end
shading('flat');axis('square'); axis('tight');
xlabel('Latitude','fontsize',16);
ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse')

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);
    
