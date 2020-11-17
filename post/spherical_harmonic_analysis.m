% Analysis spherical harmonic data
clear;
test_case  = 'drake';
run_id     = '2layer_J6';
type       = 'baroclinic_2';
%type       = 'barotropic';
cp_id      = '0030';
local      = false;
machine    = 'if.mcmaster.ca';
%machine    = 'mac';
%% Plot data
file_base   = [run_id '_' cp_id '_' type];
remote_file = ['~/hydro/' test_case '/' file_base];
local_file  = ['~/hydro/' test_case '/' file_base];
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
if ~strcmp(machine,'mac') 
    unix (sprintf(scp_cmd));
end

fid = fopen(local_file);
data = fread(fid,'single');
N = round(-3/2 + sqrt(2*(numel(data)+1)));
data = reshape (data, N+1, N/2+1)';
dmin = min(min(data));
dmax = max(max(data));
fprintf('Minimum value of data = %8.4e\n', dmin);
fprintf('Maximum value of data = %8.4e\n', dmax);

lon  = (-N/2:N/2) * 360/N;
lat  = (-N/4:N/4) * 360/N;

% Target area for local spectral analysis
if local
%     % Vortical region at southern edge of land mass
%     lat0   = -40;
%     lon0   =  15;
%     theta0 =  10;
    
    % Vortical region at 45 N
    lat0   =  45;
    lon0   =  20;
    theta0 =  20;
    
%     % Laminar region
%     lat0   =  0;
%     lon0   = -100;
%     theta0 =  30;
    
    lat_min = lat0 - theta0;
    lat_max = lat0 + theta0;
    lon_min = lon0 - theta0;
    lon_max = lon0 + theta0;
else
    lat_min = min(lat);
    lat_max = max(lat);
    lon_min = min(lon);
    lon_max = max(lon);
end
ax = [lon_min lon_max lat_min lat_max];

c_scale = linspace(dmin, dmax, 100); 
v_title = 'vorticity';
smooth  = 0;
lines   = 0;
shift   = 0;
plot_lon_lat_data(data, lon, lat, c_scale, v_title, smooth, shift, lines);
axis(ax);
%% Plot spectrum
%figure;
file_base   = [run_id '_' cp_id '_' type '_spec'];
remote_file = ['~/hydro/' test_case '/' file_base];
local_file  = ['~/hydro/' test_case '/' file_base];
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
if ~strcmp(machine,'mac') 
    unix (sprintf(scp_cmd));
end
pspec = load(local_file);
pspec(:,2) = pspec(:,2)./pspec(:,1).^2; % convert vorticity spectrum to energy spectrum integrated over shells
loglog(pspec(:,1),pspec(:,2),'b-','linewidth',1.5,'DisplayName','global');hold on;grid on;
p = -3; 
nspec=numel(pspec(:,1));
k1 = 10; k2 = round(0.3*nspec); knorm = 100;
loglog(pspec(k1:k2,1),0.9*pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'m-','linewidth',1.5,'DisplayName','k^{-3}');

% p = -8;
% nspec=numel(pspec(:,1));
% k1 = 250; k2 = round(nspec);knorm=0.8*nspec;
% loglog(pspec(k1:k2,1),1.0*pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'r-','linewidth',1.5,'DisplayName','k^{-8}');

if local
    file_base   = [run_id '_' cp_id '_' type '_local_spec'];
    remote_file = ['~/hydro/' test_case '/' file_base];
    local_file  = ['~/hydro/' test_case '/' file_base];
    scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
    if ~strcmp(machine,'mac')
        unix (sprintf(scp_cmd));
    end
    lspec = load(local_file);
    lspec(:,2) = lspec(:,2)./lspec(:,1).^2;
    loglog(lspec(:,1),lspec(:,2),'g-','linewidth',1.5,'DisplayName','local');hold on; grid on;
    
    nspec=numel(lspec(:,1));
%     p = -2.5;
%     nspec=numel(lspec(:,1));
%     k1 = 5; k2 = nspec;%round(0.2*nspec); 
%     knorm=50;
%     loglog(lspec(k1:k2,1),lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'k-','linewidth',1.5,'DisplayName','k^{-2.5}');
    p = -4;
    k1 = 5; k2 = round(nspec); 
    knorm=50;
    loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r-','linewidth',1.5,'DisplayName','k^{-4}');
end
%


titl = [type ' power spectrum'];
xlabel('l');ylabel('power');title(titl);legend;
set (gca,'fontsize',18);
