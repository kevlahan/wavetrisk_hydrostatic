function plot_lon_lat_data(s, lon, lat, c_scale, v_title, smooth, shift, lines)
% Plot longitude-latitude data
% Input:
% s       = data
% lon     = longitude coordinates
% lat     = latitude coordinates
% c_scale = colour scale for data
% v_title = name for data
% shift   = 0: left boundary at -180 longitude, 1: left boundary at zero longitude

v_yticks = [-90 -60 -30 0 30 60 90];
if shift % shift left boundary to zero longitude
    s = circshift(s,round(size(s,2)/2),2);
    axis([0 360 -90 90]);
    v_xticks = [0 30 60 90 120 150 180 210 240 270 300 330 360];
    lon = lon+180;
else
    axis([-180 180 -90 90]);
    v_xticks = [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180];
end

if smooth
    s = smooth2a(s,2,2);
end
[~,h]=contourf(lon,lat,s,100);
if not(lines)
    set(h,'LineColor','none')
end
J = customcolormap_preset('red-white-blue');
colormap(J)
%colormap(J(numel(c_scale)-1));
clim([min(c_scale) max(c_scale)]);c=colorbar;

c.Label.String=v_title; c.Label.FontSize=18;%c.YTick=c_scale;
axis('equal'); 
set(gca,'ytick',v_yticks,'yticklabels',{'90S','60S','30S','0','30N','60N','90N'});
if shift
    %set(gca,'xtick',v_xticks,'xticklabels',{'0', ' ', '60E', ' ', '120E', ' ', '180', ' ', '120W', ' ', '60W', ' '});
    set(gca,'xtick',v_xticks,'xticklabels',{'0', '30 ', '60', '90 ', '120', '150 ', '180', ...
        '210 ', '240', '270 ', '300', '330', '360'});
else
    set(gca,'xtick',v_xticks,'xticklabels',{'180W ', '150W ', '120W', '90W ', '60W', '30W ', '0', ...
        '30E ', '60E', '90E', '120E', '150E', '180E'});
end
xlabel('Longitude');
ylabel('Latitude' );
set(gca,'FontSize',19);








