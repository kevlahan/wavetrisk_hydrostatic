itype     = 'max_wave_height';  % field to plot
%itype     = 'arrival_time';  % field to plot
iwrite     = 12;     % time to plot
smooth    = false;  % smooth data over two points in each direction
shift     = true;   % shift left boundary to zero longitude
lines     = false;   % plot lines

test_case = 'tsunami'; run_id = 'tsunami'; run_dir = '';
file_base = [run_id '.' int2str(60000+iwrite)];
pathid = ['/Users/kevlahan/hydro/' test_case '/' run_dir];
file_tar = ['tar ' 'xf ' pathid file_base '.tgz'];
disp(['Uncompressing file ' pathid file_base '.tgz']);
system(file_tar);

% Load coordinates
lon = load([file_base '20']);
lat = load([file_base '21']);

ax = [0 360 -90 90];

if (strcmp(itype,'max_wave_height')) % Plot surface pressure data
    s_ll = load([file_base '11']);
    c_scale = -0.01:20:0.01;
    v_title = 'Maximum wave height (m)';
elseif (strcmp(itype,'arrival_time')) % Plot surface pressure data
    s_ll = load([file_base '12']);
    c_scale =  0:20:300;
    v_title = 'First arrival time (minutes)';
end

fprintf('Minimum value of variable %s = %8.4e\n', itype, min(min(s_ll)));
fprintf('Maximum value of variable %s = %8.4e\n', itype, max(max(s_ll)));
figure
plot_lon_lat_data(s_ll, lon, lat, c_scale, v_title, smooth, shift, lines)
axis(ax)