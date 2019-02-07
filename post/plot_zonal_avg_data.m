function plot_zonal_avg_data(s, lat, P_z, c_scale, v_title, smooth, lines, unif_grid)
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
    if (smooth)
        P_unif = smooth2a(P_unif,2,2);
    end
    [C,h]=contourf(lat_unif, P_unif, griddata(lat,P_z,s,lat_unif,P_unif), c_scale);
else
    if (smooth)
        P_z = smooth2a(P_z,2,2);
    end
    [C,h]=contourf(lat,P_z,s, c_scale);
end

if not(lines)
    set(h,'LineColor','none')
end

%colormap(jet(numel(c_scale)-1));
set(0,'defaulttextinterpreter','latex')

colormap(jet);
caxis([min(c_scale) max(c_scale)]);
colormap(jet(numel(c_scale)-1));caxis([min(c_scale) max(c_scale)]);c=colorbar;
c.Label.String=v_title;c.Label.FontSize=16;c.YTick=c_scale;
shading('flat');axis('square'); axis('tight');yticks(v_yticks);
xlabel('Latitude','fontsize',16);ylabel('$P/P_S$','fontsize',16);
set(gca,'xtick',v_xticks,'xticklabels',{'90S','60S','30S','0','30N','60N','90N'});
set(gca,'Ydir','reverse');set(gca,'FontSize',16);

