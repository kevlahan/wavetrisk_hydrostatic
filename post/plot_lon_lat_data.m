% Plot longitude-latitude data
itime = '001';

% Extract files
file_base = 'fort.3';
file_tar = ['tar ' 'xvf ' '/Users/kevlahan/hydro/' file_base itime '.tgz'];
system(file_tar);

% Load coordinates
lon = load([file_base itime '21']);
lat = load([file_base itime '22']);

% Load mass
s = load([file_base itime '01']);

% Plot data
pcolor(lon,lat,s);
c=colorbar;c.Label.String='Mass density';c.Label.FontSize=12;
shading('flat');axis('equal');axis([0 360 -90 90]); axis('tight'); 
xlabel('Longitude','fontsize',16);
ylabel('Latitude','fontsize',16);

% Erase extracted files
file_erase = ['\rm ' file_base '*'];
system(file_erase);
    
