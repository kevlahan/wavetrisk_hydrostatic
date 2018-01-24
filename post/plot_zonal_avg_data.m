function plot_zonal_avg_data(s, lat, P_z, c_scale, v_title, unif_grid)
% Plot zonally averaged data
%
% Input:
% s       = data
% lat     = latitude coordinates
% P_z     = pressure coordinates
% c_scale = colour scale for data
% v_title = name for data

v_xticks = [-90 -60 -30 0 30 60];
v_yticks = [0 0.25 0.5 0.75 1];
figure;
if (unif_grid) %Interpolate onto uniform pressure grid
    P_unif = linspace(min(P_z),max(P_z),numel(P_z));
    [lat_unif,P_unif] = meshgrid(lat,P_unif);
    contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), c_scale);
else
    contourf(lat,P_z,s);
end
colormap(jet(numel(c_scale)-1));caxis([min(c_scale) max(c_scale)]);c=colorbar;c.Label.String=v_title;c.Label.FontSize=16;
shading('flat');axis('square'); axis('tight');xticks(v_xticks);yticks(v_yticks);
xlabel('Latitude','fontsize',16);ylabel('P/P_S','fontsize',16);
set(gca,'Ydir','reverse');set(gca,'FontSize',16);