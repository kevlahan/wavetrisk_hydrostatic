% Analyse spherical harmonic data
clear;
test_case  = 'drake';
%run_id     = '2layer_slow'; cp_id = '0015';
%run_id     = '1layer_J6'; cp_id = '0015';
%run_id     = '2layer_J6'; cp_id = '0015';
%run_id     = '1layer_fast'; cp_id = '0015';
run_id     = '2layer_fast'; cp_id = '0015';
%type       = 'baroclinic_2';
type       = 'total_2';
%type       = 'barotropic';
local      = false;
machine    = 'if.mcmaster.ca';
%machine    = 'mac';
new_transfer = true;

% Parameters
scale_earth =  6;
if strcmp(run_id,'1layer_fast')
    scale_omega =  1;
    uwbc        =  1.0;
elseif strcmp(run_id,'2layer_fast')
    scale_omega =  1;
    uwbc        =  1.0;
elseif strcmp(run_id,'1layer_J6')
    scale_omega =  6;
    uwbc        =  1.8;
elseif strcmp(run_id,'2layer_J6')
    scale_omega =  6;
    uwbc        =  1.8;
elseif strcmp(run_id,'2layer_slow')
    scale_omega =  24;
    uwbc        =  1.8;
end

theta       =  45; % latitude at which to calculate f0 and beta
omega       =  7.29211e-5/scale_omega;
radius      =  6371.229e3/scale_earth;   
g           =  9.80616;
drho        = -4;
ref_density =  1028;
H1          =  3e3;
H2          =  1e3;
H           =  H1 + H2;
f0          =  2*omega*sin(deg2rad(theta));
beta        =  2*omega*cos(deg2rad(theta))/radius;
c0          =  sqrt(g*H);
c1          =  sqrt (g*abs(drho)/2/ref_density * H2*(H-H2)/H);
visc        =  8.25e-1;

% Lengthscales (km)
lambda0    = c0/f0;             % external radius of deformation
lambda1    = c1/f0;             % internal radius of deformation
deltaSM    = uwbc/f0;           % submesoscale
deltaI     = sqrt(uwbc/beta);   % inertial layer
deltaM     = (visc/beta)^(1/3); % Munk layer 
%% Plot data
file_base   = [run_id '_' cp_id '_' type];
remote_file = ['~/hydro/' test_case '/' file_base];
local_file  = ['~/hydro/' test_case '/' file_base];
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
if ~strcmp(machine,'mac') && new_transfer 
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
    % Vortical region at southern edge of land mass
    lat0   = -40;
    lon0   =  20;
    theta0 =  30;

%     % Vortical region at southern edge of land mass 2layer_fill
%      lat0   = -40;
%      lon0   =  10;
%      theta0 =  50;
    
%     % Vortical region at equator
%     lat0   = 0;
%     lon0   = 20;
%     theta0 = 15;
    
%     % Vortical region at 45 N
%     lat0   =  45;
%     lon0   =  20;
%     theta0 =  20;
    
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
scales = 2*pi*radius/1e3./sqrt(pspec(:,1).*(pspec(:,1)+1)); % equivalent length scale (Jeans relation)

if strcmp(type,'barotropic')
    % Power spectrum
    loglog(scales(:,1),pspec(:,2),'b-','linewidth',3,'DisplayName',type);hold on;grid on;
    
    % Power laws
    if strcmp(run_id,'2layer_slow')
        powerlaw (scales, pspec(:,2), 2000, 20, -2.7, 'm--')
        powerlaw (scales, pspec(:,2),   20,  2, -6, 'c--')
    elseif strcmp(run_id,'1layer_fast')
        powerlaw (scales, pspec(:,2), 700, 20, -3, 'm--')
        powerlaw (scales, pspec(:,2),  20,  2, -5, 'c--')
    elseif strcmp(run_id,'2layer_fast')
        powerlaw (scales, pspec(:,2), 700, 20, -3, 'm--')
        powerlaw (scales, pspec(:,2),  20,  2, -5, 'c--')
    elseif strcmp(run_id,'1layer_J6')
        powerlaw (scales, pspec(:,2), 700, 20, -3, 'm--')
        powerlaw (scales, pspec(:,2),  20,  2, -6, 'c--')
    elseif strcmp(run_id,'2layer_J6')
        powerlaw (scales, pspec(:,2), 700, 20, -3, 'm--')
        powerlaw (scales, pspec(:,2),  20,  2, -6, 'c--')
    end
end

if strcmp(type,'baroclinic_1') || strcmp(type,'baroclinic_2')
    % Power spectrum
    loglog(scales(:,1),pspec(:,2),'r-','linewidth',3,'DisplayName',type);hold on;grid on;
    
    % Power laws
    if strcmp(run_id,'2layer_slow')
        powerlaw (scales, pspec(:,2), 1500, 15, -3, 'b--')
    elseif strcmp(run_id,'2layer_fast')
        powerlaw (scales, pspec(:,2), 2000, 10, -5/3, 'b--')
    elseif strcmp(run_id,'2layer_J6')
        powerlaw (scales, pspec(:,2), 2000, 10, -5/3, 'b--')
        powerlaw (scales, pspec(:,2),  20,  2, -6, 'c--')
    end
end

if strcmp(type,'total_1') || strcmp(type,'total_2')
    % Power spectrum
    loglog(scales(:,1),pspec(:,2),'r-','linewidth',3,'DisplayName',type);hold on;grid on;
    
    % Power laws
    if strcmp(run_id,'2layer_slow')
        powerlaw (scales, pspec(:,2), 1500, 15, -3, 'b--')
    elseif strcmp(run_id,'2layer_fast')
        powerlaw (scales, pspec(:,2), 2000, 10, -3, 'b--')
        powerlaw (scales, pspec(:,2),  20,  2, -5, 'c--')
    elseif strcmp(run_id,'2layer_J6')
        powerlaw (scales, pspec(:,2), 2000, 10, -5/3, 'b--')
        powerlaw (scales, pspec(:,2),  20,  2, -6, 'c--')
    end
end

% Local power spectra
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
    lscales = 2*pi*radius/1e3./sqrt(pspec(:,1).*(pspec(:,1)+1)); % equivalent length scale (Jeans relation)

    loglog(lscales,lspec(:,2),'g-','linewidth',2,'DisplayName','local');hold on; grid on;
    
    nspec=numel(lspec(:,1));
end

xlabel("l (km)");ylabel("S(l)");

if strcmp(run_id,'2layer_slow')
    title("\Omega = \Omega_{Earth}/24");
elseif strcmp(run_id,'1layer_fast') || strcmp(run_id,'2layer_fast')
    title("\Omega = \Omega_{Earth}");
elseif strcmp(run_id,'1layer_J6')
    title("\Omega = \Omega_{Earth}/6");
elseif strcmp(run_id,'2layer_J6')
    title("\Omega = \Omega_{Earth}/6");
end
set (gca,'Xdir','reverse');legend;
if strcmp(type,'barotropic')
    axis([1 1e4 1e-12 1e1]);
elseif strcmp(type,'baroclinic_1') || strcmp(type,'baroclinic_2')
    axis([1 1e4 1e-13 1e-0]);
end
set (gca,'fontsize',20);
plot_scale(lambda0/1e3,"\lambda_0");
plot_scale(lambda1/1e3,"\lambda_1");
plot_scale(deltaSM/1e3,"\delta_{SM}");
%plot_scale(deltaM/1e3,"\delta_{M}");

%%
function powerlaw (scales, power, s1, s2, p, col)
% Plots power law -p between scales s1 and s2 (s2 > s1) with color col
[~,k1] = min(abs(scales-s1)); [~,k2] = min(abs(scales-s2)); knorm = round((k1+k2)/2);
str = "k^{" + compose("%1.2f",p) + "}";
loglog(scales(k1:k2),scales(k1:k2).^(-p) * power(knorm)/scales(knorm)^(-p),col,'linewidth',3,'DisplayName',str);
end

function plot_scale (scale,name)
y = ylim;
x = scale;
h=loglog([x x], y, 'k','linewidth',1.5,'linewidth',2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(0.85*scale,10*y(1),name,"fontsize",16)
end
