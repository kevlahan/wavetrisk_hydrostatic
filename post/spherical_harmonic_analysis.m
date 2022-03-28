%% Analyse spherical harmonic data
test_case  = 'drake';
type       = 'barotropic';
local      = true;
machine    = 'if.mcmaster.ca';
%machine    = 'niagara.computecanada.ca';

% Parameters
scale_earth =  6;
scale_omega =  6;
uwbc        =  1.0;

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

file_base   = [run_id '_' cp_id '_' type];
remote_file = ['~/hydro/' test_case '/' file_base];
local_file  = ['~/hydro/' test_case '/' file_base];
%% Load local region data
if ~strcmp(machine,'mac') 
    scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
    unix (sprintf(scp_cmd));
end
%% Plot local region

local      = true;
region     = 'mid';
run_id     = "2layer_normal";
cp_id      = "0031";
type       = "barotropic";
local_file = run_id + "_" + cp_id + "_" + type;

fid = fopen(local_file);
data = fread(fid,'double'); 
N = round(-3/2 + sqrt(2*(numel(data)+1)));
%data = reshape (data, N+1, N/2+1)';
data = reshape (data, N+1, N/2+1)';
dmin = min(min(data));
dmax = max(max(data));
fprintf('Minimum value of data = %8.4e\n', dmin);
fprintf('Maximum value of data = %8.4e\n', dmax);

lon  = (-N/2:N/2) * 360/N;
lat  = (-N/4:N/4) * 360/N;

% Target area for local spectral analysis

if local
    if strcmp(region,'southern') % Vortical region at southern edge of land mass
        lat0   = -40;
        lon0   =  20;
        theta0 =  30;
    elseif strcmp(region,'equator') % Vortical region at equator
        lat0   = 0;
        lon0   = 35;
        theta0 = 20;
    elseif strcmp(region,'mid') % Vortical region at 45 N
        lat0   =  45;
        lon0   =  35;
        theta0 =  20;
    elseif strcmp(region,'laminar') % Laminar region
        lat0   =  0;
        lon0   = -100;
        theta0 =  20;
    end
    
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

vort_limit = 3e-5;

% Color vorticity beyond vort_limit to maximum colors
data(data > vort_limit)  =  vort_limit;
data(data < -vort_limit) = -vort_limit;

c_scale = linspace(-vort_limit, vort_limit, 100); 

if strcmp(type,"barotropic")
    v_title = 'Barotropic vorticity';
elseif strcmp(type,"baroclinic")
    v_title = 'Baroclinic vorticity';
end
smooth  = 0;
lines   = 0;
shift   = 0;
plot_lon_lat_data(data(1:4:end,1:4:end), lon(1:4:end), lat(1:4:end), c_scale, v_title, smooth, shift, lines);
axis(ax);
%% Load spectrum data
% file_base   = [run_id '_' cp_id '_' type '_spec'];
% remote_file = ['~/hydro/' test_case '/' file_base];
% local_file  = ['~/hydro/' test_case '/' file_base];
% scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
% if ~strcmp(machine,'mac') 
%     unix (sprintf(scp_cmd));
% end

% Plot spectrum
run_id      = "2layer_slow";
scale_omega = 6;
uwbc        = 1.5; 
local       = false;
average     = true;

cp_id       = "0014";
type        = "barotropic";
%type        = "baroclinic_2";
%type        = "total_2";

p1           = -5/3; % comparison power spectrum
range1       = [1000 10];    
range1       = [200 5];    
%range1       = [50 4]; 

if average
    if local
        local_file  = run_id + "_" + type + "_local_spec"; % average
    else
        local_file  = run_id + "_" + type + "_spec"; % average
    end
else
    if local
        local_file  = run_id + "_" + cp_id + "_" + type + "_local_spec";
    else
        local_file  = run_id + "_" + cp_id + "_" + type + "_spec";
    end
end

% Physical parameters of simulation
scale_earth = 6;
theta       =  45; % latitude at which to calculate f0 and beta
omega       =  7.29211e-5/scale_omega;
radius      =  6371.229e3/scale_earth; 
visc        =  3.94;
g           =  9.80616;
drho        = -4;
ref_density =  1028;
H1          =  3e3;
H2          =  1e3;
H           =  H1 + H2;
f0          =  2*omega*sin(deg2rad(theta));
beta        =  2*omega*cos(deg2rad(theta))/radius;
r_b         =  1e-7;
c0          =  sqrt(g*H);
c1          =  sqrt (g*abs(drho)/ref_density * H2*(H-H2)/H);

if strcmp(type,'barotropic')
    name_type = "Barotropic";
elseif strcmp(type,'baroclinic_1')
    name_type = "Lower layer baroclinic";
elseif strcmp(type,'baroclinic_2')
    name_type = "Top layer baroclinic";
elseif strcmp(type,'total_1')
    name_type = "Lower layer total";
elseif strcmp(type,'total_2')
    name_type = "Top layer total";
end

% Lengthscales (km)
lambda0    = c0/f0;             % external radius of deformation
lambda1    = c1/f0;             % internal radius of deformation
deltaS     = r_b / beta;        % Stommel layer
deltaSM    = uwbc/f0;           % submesoscale
deltaI     = sqrt(uwbc/beta);   % inertial layer
deltaM     = (visc/beta)^(1/3); % Munk layer 

pspec = load(local_file);
pspec(:,2) = pspec(:,2)./pspec(:,1).^2; % convert vorticity spectrum to energy spectrum integrated over shells
scales = 2*pi*radius/1e3./sqrt(pspec(:,1).*(pspec(:,1)+1)); % equivalent length scale (Jeans relation)

% Power spectrum
loglog(scales(:,1),pspec(:,2),'b-','linewidth',3,'DisplayName',name_type);hold on;grid on;

% Log-law comparison
powerlaw (scales, pspec(:,2), range1, p1, 'r--')

axis([4e0 1e4 1e-10 1e0]);
set (gca,'fontsize',20);
xlabel("\lambda (km)");ylabel("S(\lambda)");
set (gca,'Xdir','reverse');legend;
plot_scale(lambda0/1e3,"\lambda_0");
plot_scale(lambda1/1e3,"\lambda_1");
plot_scale(deltaSM/1e3,"\delta_{SM}");
%plot_scale(deltaS/1e3,"\delta_{S}");
%plot_scale(deltaM/1e3,"\delta_{M}");

%%
function powerlaw (scales, power, range, p, col)
% Plots power law -p between scales s1 and s2 (s2 > s1) with color col
[~,k1] = min(abs(scales-range(1))); [~,k2] = min(abs(scales-range(2))); knorm = round((k1+k2)/2);

if abs(p+5/3) < 1e-2
    str = "k^{"+ "-"+ 5 +"/" + 3 +"}";
elseif mod(p,1) ~= 0
    str = "k^{" + compose("%1.1f",p) + "}";
else
    str = "k^{" + compose("%1.0f",p) + "}";
end

loglog(scales(k1:k2),scales(k1:k2).^(-p) * power(knorm)/scales(knorm)^(-p),col,'linewidth',3,'DisplayName',str);
end

function plot_scale (scale,name)
y = ylim;
x = scale;
h=loglog([x x], y, 'k','linewidth',1.5);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
text(0.92*scale,10*y(1),name,"fontsize",16)
end
