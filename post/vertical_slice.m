clear
% test_case = 'seamount';
% run_id    = 'sea_drho_3_nonadapt';

test_case = 'upwelling';
run_id    = 'upwelling';

%machine   = 'if.mcmaster.ca';
machine   = 'cherry';

% Transfer data
directory   = ['~/hydro/' test_case];
file_base   = [run_id '.5'];

if ~strcmp(machine,'cherry')
remote_file = [directory '/' file_base '.tgz'];
local_file  = remote_file;
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
unix (sprintf(scp_cmd));
file_tar = ['tar ' 'xf ' local_file ' -C ' directory];
disp(['Uncompressing file ' local_file]);
system(file_tar);
end
%%
zlevels = 16;

% Load coordinates
lon = fread(fopen([directory '/' file_base '.50']),'double'); Nlon = numel(lon);
lat = fread(fopen([directory '/' file_base '.51']),'double'); Nlat = numel(lat); 

xlat = fread(fopen([directory '/' file_base '.52']),[Nlat,2],'double'); xlat = [xlat ; -xlat(1,:)];
xlon = fread(fopen([directory '/' file_base '.53']),[Nlon,2],'double');

zlat = fread(fopen([directory '/' file_base '.54']),'double'); zlat = reshape(zlat,[Nlat,zlevels,2]); 
zlon = fread(fopen([directory '/' file_base '.55']),'double'); zlon = reshape(zlon,[Nlon,zlevels,2]); 

lat_slice = fread(fopen([directory '/' file_base '.56']),'double');lat_slice = reshape(lat_slice,[Nlat,zlevels,4]);
lon_slice = fread(fopen([directory '/' file_base '.57']),'double');lon_slice = reshape(lon_slice,[Nlon,zlevels,4]);
%% Plot vertical grid
hf1=figure(1);
plot(lat,zlat(:,:,1)','k-','linewidth',1.2); axis([30 60 -150 0]);
xlabel('latitude'); ylabel('z (km)'); hold on;
plot([xlat(:,1) xlat(:,1)], ylim, 'k-','linewidth',1.2);
%set(gca,'DataAspectRatio',[10 1 1])
%pos = get(hf1,'position');
%set(hf1,'position',[pos(1:2)/4 pos(3:4)*2])
set(gca,'fontsize',18);hold off
%% Plot data
hf2=figure(2);
for i = 1:Nlat
    for k = 1:zlevels
        x = [xlat(i,1) xlat(i,2) xlat(i,2) xlat(i,1)];
        y = [zlat(i,k,1) zlat(i,k,1) zlat(i,k,2) zlat(i,k,2)];
        patch(x,y,lat_slice(i,k,3),'LineStyle','none');
    end
end
hcb=colorbar;hcb=colorbar; hcb.Label.String = "Density perturbation (kg/m^3)";
xlabel('latitude'); ylabel('z (km)');colormap(flipud(jet)) 
axis([35 55 -150 0]);hold on;
set(gca,'fontsize',18);
%set(gca,'DataAspectRatio',[10 1 1])
%pos = get(hf2,'position');
%set(hf2,'position',[pos(1:2)/4 pos(3:4)*2])
%hcb.FontSize=10;



