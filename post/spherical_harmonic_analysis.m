%% Load data
clear
test_case = 'drake';
%machine  = 'if.mcmaster.ca';
machine   = 'niagara.computecanada.ca';

%run_id = '1layer_slow' ; type = 'barotropic_spec';       name_type = '1 layer'; cp_id = '0031';
%run_id = '2layer_slow' ; type = 'total_2_spec';          name_type = '2 layer'; cp_id = '0031';
%run_id = '6layer_test' ; type = 'barotropic_curlu_spec'; name_type = '6 layer'; cp_id = '0038';
run_id = '12layer_test'; type = 'barotropic_curlu_spec'; name_type = '12 layer'; cp_id = '0008';
avg = true;

local_file = load_data(test_case, run_id, cp_id, type, machine, avg);

%% Plot energy spectra
col_spec  = 'b-';  % colour for energy spectrum
col_power = 'r--'; % colour for power law

% Physical parameters of simulation
visc        =  0.493;
uwbc        =  0.8; 
scale_omega =  6;
scale_earth =  6;
omega       =  7.29211e-5/scale_omega;
radius      =  6371.229e3/scale_earth; 
g           =  9.80616;
drho        = -4;
ref_density =  1028;
H1          =  3e3;
H2          =  1e3;
H           =  H1 + H2;
theta       =  45; % latitude at which to calculate f0 and beta
f0          =  2*omega*sin(deg2rad(theta));
beta        =  2*omega*cos(deg2rad(theta))/radius;
%r_b         =  1.3e-8; % two-layer
r_b         =  1e-7;
c0          =  sqrt(g*H);
%c1          =  sqrt (g*abs(drho)/ref_density * H2*(H-H2)/H); % two-layer
c1          =  5.56; % m/s

% Lengthscales (km)
lambda0    = c0/f0/1e3;             % external radius of deformation
lambda1    = c1/f0/1e3;             % internal radius of deformation
deltaS     = r_b / beta/1e3;        % Stommel layer
deltaSM    = uwbc/f0/1e3;           % submesoscale
deltaI     = sqrt(uwbc/beta)/1e3;   % inertial layer
deltaM     = (visc/beta)^(1/3)/1e3; % Munk layer 

pspec = load(local_file);
pspec(:,2) = pspec(:,2)./pspec(:,1).^2;                     % convert vorticity spectrum to energy spectrum integrated over shells
scales = 2*pi*radius/1e3./sqrt(pspec(:,1).*(pspec(:,1)+1)); % equivalent length scale (Jeans relation)

% Power spectrum
loglog(scales,pspec(:,2),col_spec,'linewidth',3,'DisplayName',name_type);hold on;grid on;

% Fit power law
range = [lambda1 deltaSM ]; % range to fit power law 
fit_indices = find(scales > range(2) & scales < range(1));
[P,S] = polyfit(log10(scales(fit_indices)),log10(pspec(fit_indices,2)),1);

st_err = sqrt(diag(inv(S.R)*inv(S.R'))*S.normr^2/S.df); % error in coefficients from covariance matrix of P

fprintf('\nFitted power law is %.2f +/- %.2f\n',-P(1), st_err(1));

% Plot fit
powerlaw (scales, 1.5*pspec(:,2), range, -P(1), col_power)

axis([4e0 1e4 1e-10 1e0]);
set (gca,'fontsize',20);
xlabel("\lambda (km)");ylabel("S(\lambda)");
set (gca,'Xdir','reverse');legend;
%plot_scale(lambda0,"\lambda_0");
plot_scale(lambda1,"\lambda_1");
plot_scale(deltaSM,"\delta_{SM}");
%plot_scale(deltaI,"\delta_{I}");
%plot_scale(deltaS,"\delta_{S}");
%plot_scale(deltaM,"\delta_{M}");
%% Plot local region
local  = true;
region = 'mid';

fid = fopen(local_file);
data = fread(fid,'double'); 
N = round(-3/2 + sqrt(2*(numel(data)+1)));
data = reshape (data, N+1, N/2+1)';
dmin = min(min(data)); dmax = max(max(data));
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

function [local_file] =  load_data (test_case, run_id, cp_id, type, machine, avg)
if avg % average spectrum
    file_base = [run_id '_' type];
else
    file_base = [run_id '_' cp_id '_' type];
end
remote_file = ['~/hydro/' test_case '/' file_base];
local_file  = ['~/hydro/' test_case '/' file_base];
scp_cmd     = ['scp ' machine ':' remote_file ' ' local_file];
if ~strcmp(machine,'mac')
    unix (sprintf(scp_cmd));
end
end


