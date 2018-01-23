% Plot longitude-latitude data
clear all; close all;
itime      = '006';
unif_grid  = false;
n_contours = 11;

v_xticks = [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150];
v_yticks = [-90 -60 -30 0 30 60 90];

% Extract files
file_base = 'fort.3';
file_tar = ['tar ' 'xvf ' '/Users/kevlahan/hydro/' file_base itime '.tgz'];
system(file_tar);

% Load coordinates
lon = load([file_base itime '20']);
lat = load([file_base itime '21']);
P_z = load([file_base itime '22']); % Pressure-based vertical coordinates

% Plot mass density data
s = load([file_base itime '01']); 
figure(1);colormap(jet(n_contours))
contourf(lon,lat,s,n_contours);hold on; 
c=colorbar;c.Label.String='Mass density';c.Label.FontSize=12;
shading('flat');axis('equal');axis([-180 180 -90 90]); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Longitude','fontsize',16);ylabel('Latitude','fontsize',16);

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
shading('flat');axis('square'); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Latitude','fontsize',16);ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse')

% Plot temperature data
v = 270:3:303;
s = load([file_base itime '02']); 
figure(3);
contourf(lon,lat,s,n_contours);hold on; 
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Temperature (K)';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Longitude','fontsize',16);ylabel('Latitude','fontsize',16);

% Plot zonally averaged temperature
s = load([file_base itime '12']);
figure(4); 
if (unif_grid)
    %Interpolate onto uniform pressure grid
    P_unif = linspace(min(P_z),max(P_z),numel(P_z));
    [lat_unif,P_unif] = meshgrid(lat,P_unif);
    contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), v);
    else
    contourf(lat,P_z,s);
end
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Temperature (K)';c.Label.FontSize=12;
shading('flat');axis('square'); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Latitude','fontsize',16);ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse')

% Plot zonal velocity data
v = -15:5:50;
s = load([file_base itime '03']); 
figure(5);
contourf(lon,lat,s,v);hold on; 
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Zonal velocity (m/s)';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Longitude','fontsize',16);ylabel('Latitude','fontsize',16);

% Plot zonally averaged zonal velocity
s = load([file_base itime '13']);
figure(6); colormap(jet(numel(v)-1))
if (unif_grid)
    %Interpolate onto uniform pressure grid
    P_unif = linspace(min(P_z),max(P_z),numel(P_z));
    [lat_unif,P_unif] = meshgrid(lat,P_unif);
    contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), v);
else
    contourf(lat,P_z,s);
end
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Zonal velocity (m/s)';c.Label.FontSize=12;
shading('flat');axis('square'); axis('tight');
xlabel('Latitude','fontsize',16);ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse')

% Plot meridional velocity data
v = -35:5:20;
s = load([file_base itime '04']); 
figure(7);
contourf(lon,lat,s,v); 
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Meridional velocity (m/s)';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Longitude','fontsize',16);ylabel('Latitude','fontsize',16);

% Plot zonally averaged meridional velocity
s = load([file_base itime '14']);
figure(8); 
if (unif_grid)
    %Interpolate onto uniform pressure grid
    P_unif = linspace(min(P_z),max(P_z),numel(P_z));
    [lat_unif,P_unif] = meshgrid(lat,P_unif);
    contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), v);
else
    contourf(lat,P_z,s,v);
end
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Meridional velocity (m/s)';c.Label.FontSize=12;
shading('flat');axis('square'); axis('tight');
xlabel('Latitude','fontsize',16);ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse')

% Plot geopotential
v = 2400:100:3400;
s = load([file_base itime '05']); 
figure(9);
contourf(lon,lat,s,v); 
colormap(jet(numel(v)-1));caxis([min(v) max(v)]);c=colorbar;c.Label.String='Geopotential (m)';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight');xticks(v_xticks);yticks(v_yticks); 
xlabel('Longitude','fontsize',16); ylabel('Latitude','fontsize',16);

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);
    
