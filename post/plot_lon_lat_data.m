function plot_lon_lat_data(s, lon, lat, c_scale, v_title, smooth, shift)
% Plot longitude-latitude data
% Input:
% s       = data
% lon     = longitude coordinates
% lat     = latitude coordinates
% c_scale = colour scale for data
% v_title = name for data
% shift   = 0: left boundary at -180 longitude, 1: left boundary at zero longitude

axis([-180 180 -90 90]);
v_xticks = [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150];
v_yticks = [-90 -60 -30 0 30 60 90];
if shift % Shift left boundary to zero longitude
    s = circshift(s,round(size(s,2)/2),2);
    axis([0 360 -90 90]);
    v_xticks = [0 30 60 90 120 150 180 210 240 270 300 330];
    lon = lon+180;
end
if (smooth)
    s = smooth2a(s,5,5);
end
contourf(lon,lat,s,c_scale);

colormap(jet(numel(c_scale)-1));caxis([min(c_scale) max(c_scale)]);c=colorbar;c.Label.String=v_title;c.Label.FontSize=16;
axis('equal'); 
set(gca,'ytick',v_yticks,'yticklabels',{'90S','60S','30S','0','30N','60N','90N'});
if (shift)
    %set(gca,'xtick',v_xticks,'xticklabels',{'0', ' ', '60E', ' ', '120E', ' ', '180', ' ', '120W', ' ', '60W', ' '});
    set(gca,'xtick',v_xticks,'xticklabels',{'0', '30 ', '60', '90 ', '120', '150 ', '180', ...
        '210 ', '240', '270 ', '300', '330'});
end
xlabel('Longitude','fontsize',16);ylabel('Latitude','fontsize',16);set(gca,'FontSize',16);








