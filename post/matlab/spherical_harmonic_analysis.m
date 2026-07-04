%% Load data
clear; clc

machine    = "gra-dtn1.alliancecan.ca";

dir_remote = "~/proj/drake/slow_rotation/J8Z60/spectra"; dir_local = "~/hydro/drake";

% Transfer all spectrum files at once
scp_cmd = "scp "+machine+":"""+dir_remote+"/*spec"" "+dir_local;
if ~strcmp(machine,"mac")
    unix (sprintf(scp_cmd));
end
%%
Z_tanh = [ ...
    -3.8000e+03  -3.4300e+03  -3.1000e+03  -2.8000e+03  -2.5400e+03  -2.3000e+03  -2.0900e+03  -1.9000e+03  -1.7300e+03  -1.5800e+03  -1.4400e+03  ...
    -1.3200e+03  -1.2100e+03  -1.1100e+03  -1.0200e+03  -9.4600e+02  -8.7400e+02  -8.0900e+02  -7.5100e+02  -6.9800e+02  -6.5000e+02  -6.0600e+02  ...
    -5.6600e+02  -5.3000e+02  -4.9700e+02  -4.6600e+02  -4.3800e+02  -4.1200e+02  -3.8900e+02  -3.6600e+02  -3.4600e+02  -3.2600e+02  -3.0800e+02  ...
    -2.9100e+02  -2.7500e+02  -2.6000e+02  -2.4500e+02  -2.3100e+02  -2.1800e+02  -2.0500e+02  -1.9300e+02  -1.8100e+02  -1.7000e+02  -1.5900e+02  ...
    -1.4800e+02  -1.3800e+02  -1.2700e+02  -1.1700e+02  -1.0700e+02  -9.7400e+01  -8.7800e+01  -7.8300e+01  -6.8800e+01  -5.9500e+01  -5.0200e+01  ...
    -4.1000e+01 -3.1900e+01  -2.2700e+01  -1.3600e+01  -4.5400e+00];

Z_linear = [ ...
    -3.7900e+03  -3.3900e+03  -3.0400e+03  -2.7200e+03  -2.4400e+03  -2.1900e+03  -1.9600e+03  -1.7600e+03  -1.5800e+03  -1.4200e+03  -1.2800e+03  ...
    -1.1600e+03  -1.0400e+03  -9.4300e+02  -8.5300e+02  -7.7200e+02  -7.0000e+02  -6.3600e+02  -5.7800e+02  -5.2600e+02  -4.8000e+02  -4.3800e+02  ...
    -4.0100e+02  -3.6700e+02  -3.3700e+02  -3.0900e+02  -2.8500e+02  -2.6200e+02  -2.4200e+02  -2.2300e+02  -2.0700e+02  -1.9100e+02  -1.7700e+02  ...
    -1.6500e+02  -1.5300e+02  -1.4200e+02  -1.3200e+02  -1.2300e+02  -1.1400e+02  -1.0600e+02  -9.8400e+01  -9.1300e+01  -8.4600e+01  -7.8300e+01  ...
    -7.2300e+01  -6.6500e+01  -6.1000e+01  -5.5800e+01  -5.0700e+01  -4.5800e+01  -4.1000e+01  -3.6400e+01  -3.1900e+01  -2.7500e+01  -2.3100e+01  ...
    -1.8800e+01  -1.4600e+01  -1.0400e+01  -6.2200e+00  -2.0700e+00];


%% Analyze spectrum data
clc; global KM; global tanh_strat; format short e

dir           = "~/hydro/drake/J8Z60_linear/";
tanh_strat    = false;        % tanh profile (or linear) for drake case
test_case     = "drake";        % prefix for test case
level         = 8;           % resolution level
zlevels       = 60;          % number of vertical layers
type          = "flux";      % u, curlu, divu or topo
avg           = true;       % analyze average spectrum or individual spectra
power_law_fit = false;       % plot power law fit
plot_spec     = true;        % plot spectrum

if tanh_strat
    Z = Z_tanh;
else
    Z = Z_linear;
end

if tanh_strat
    layers = [1 23 54 60]; % tanh
else
    layers = [1 14 46 60]; % linear
end
layers = 1:2:zlevels;
if avg
    cp_min = 1; cp_max=cp_min; 
else
    cp_min = 0; cp_max = 0;
end

% Set physical parameters
if strcmp(test_case,"drake")
    [H, lambda0, lambda1, deltaS, deltaSM, deltaRh, deltaM, radius, Z] = params();
else
    radius = 6371e3;
end

% Power law fitting ranges
range = zeros(zlevels,2);
if tanh_strat
    range( 1:15,1) = 250; range(1 :15,2) = 80;
    range(16:20,1) = 250; range(16:20,2) = 80;
    range(21:35,1) = 260; range(21:35,2) = 110;
    range(25:40,1) = 180; range(25:60,2) = 100;
    range(41:49,1) = 180; range(41:49,2) = 90;
    range(50:60,1) = 160; range(50:60,2) = 60;
else
    range( 1:30,1) = 170; range( 1:30,2) = 50;
    range(31:60,1) = 140; range(31:60,2) = 50;
end
if strcmp(type,"topo")
    range(:,1) = 1000; range(:,2) = 100;
end

plot_scales = false ;     % plot length scales
col_spec    = "b-";      % colour for energy spectrum
col_power   = "k-";      % colour for power law

run_id = test_case+"J"+num2str(level,'%1.1d')+"Z"+num2str(zlevels,'%2.2d');  % naming convention
tmpdir = dir+"temp_spec"; mkdir(tmpdir)

% Ensure each plot uses different colors
figure('Visible','on','Units','pixels','Position',[100 700 800 600]);pbaspect([4 3 1]);
%figure('Visible','on','Units','pixels','Position',[100 700 1600 1200]);pbaspect([4 3 1]);

ncurve = numel(layers);
if ncurve <= 8
    colors = [                 % Okabe-Ito color blind palette
        0.0000 0.4470 0.6980   % blue
        0.8359 0.3686 0.0000   % vermillion
        0.0000 0.6196 0.4510   % bluish green
        0.9412 0.8941 0.2588   % yellow
        0.3373 0.7059 0.9137   % sky blue
        0.8000 0.4745 0.6549   % reddish purple
        0.9020 0.6235 0.0000   % orange
        0.0000 0.0000 0.0000   % black
        ];
else
    colors = parula(ncurve);
end
ax = gca;
ax.XScale = 'log';
if any(type == ["flux", "transfer"])
    ax.YScale = 'lin';
else
    ax.YScale = 'log';
end

hp = gobjects(ncurve,1);
icurve = 0;

if power_law_fit
    pow_law = zeros(cp_max-cp_min+1,zlevels);
end

ymin = 1e16; ymax = -1e16; 
cp_idx = 0;
for cp_id = cp_min:cp_max
    cp_idx = cp_idx + 1;
    if ~avg
        cp = compose("%04d",cp_id);
        tgzfile = dir+run_id+"_spec.tgz"; 
        if power_law_fit
            fprintf('\nPower law exponents for checkpoint %d\n', cp_id)
        end
    else
        tgzfile = dir+run_id+"_spec.tgz";
    end 
    untar(tgzfile, tmpdir); % uncompress spectrum tar file

    if power_law_fit
        fprintf("Layer     p")
    end

    for zlev = fliplr(layers)

        % Load data
        k = compose("%04d",zlev);
        if strcmp(type,"topo")
            file_base = run_id+"_"+cp+"_"+type;
        else
            if avg % average spectrum
                file_base = run_id+"_"+k+"_"+type;
            else
                file_base = run_id+"_"+cp+"_"+k+"_"+type;
            end
        end
        pspec = load (tmpdir+"/"+file_base+"_spec", '-ascii');

        % Skip mean mode
        scales = 2*pi*radius/1e3./sqrt(pspec(2:end,1).*(pspec(2:end,1)+1)); % equivalent length scale (Jeans relation)
        power  = pspec(2:end,2); % power (variance) spectrum (equivalent to usual amplitude squared integrated over shells)
        rms    = pspec(2:end,3); % RMS power spectrum used in geodesy rms = sqrt(power/(2l+1))

        disp([scales(20) zlev Z_tanh(zlev) sign(power(20))])

        % Plot energy spectra

        % Convert vorticity/divergence spectrum to equivalent velocity
        % spectrum
        if any(strcmp(type, ["divu","curlu"]))
            power = power./pspec(2:end,1).^2;
        end
           
        % Fit power law
        if power_law_fit
            fit_indices = find(scales > range(zlev,2) & scales < range(zlev,1));
            [P,S] = polyfit(log10(scales(fit_indices)), log10(power(fit_indices)), 1);
            st_err = sqrt(diag(inv(S.R)*inv(S.R'))*S.normr^2/S.df); % error in coefficients from covariance matrix of P

            %fprintf("\n %3.0f    %.2f +/- %.2f", zlev, -P(1), st_err(1));
            fprintf("\n %3.0f    %.2f", zlev, -P(1)); % no fit error
            pow_law (cp_idx,zlev) = -P(1);
        end

        if plot_spec
            icurve = icurve + 1;

            if strcmp(test_case,"drake")
                if power_law_fit
                    name_type = "z = "+compose('%5.0f',Z(zlev))+" m, p = "+compose('%2.1f', -P(1));
                else
                    name_type = "z = "+compose('%5.0f',Z(zlev))+" m";
                end

                if any(type == ["flux", "transfer"])
                    y = power/max(abs(power));
                    h = semilogx(scales, y, ...
                        "LineWidth", 2, ...
                        "DisplayName", name_type, ...
                        "Color", colors(icurve,:));
                else
                    y = power;
                    h = loglog(scales, y, ...
                        "LineWidth", 2, ...
                        "DisplayName", name_type, ...
                        "Color", colors(icurve,:));
                end
            else
                name_type = "J = "+compose('%1.0d',level);

                h = loglog(scales, power, ...
                    "LineWidth", 2, ...
                    "DisplayName", name_type, ...
                    "Color", colors(icurve,:));
            end
            hold on
            grid on
            hp(icurve) = h;

            if power_law_fit % plot fit
                powerlaw (scales, 1.2*power, [range(zlev,1) range(zlev,2)], -P(1))
            end

            if strcmp(ax.YScale,'log')
                ymin = min (ymin, 10^(floor(log10(min(power)))));
                ymax = max (ymax, 10^(ceil(log10(max(power)))));
            else
                ymin = min (ymin, 1.5*min(y));
                ymax = max (ymax, 1.5*max(y));
            end
        end
    end
end
rmdir(tmpdir, 's'); % delete temporary directory

fprintf("\n")

xmin = 10^(floor(log10(min(scales))));
xmax = 10^(ceil(log10(max(scales))));
xlim([xmin xmax])
ylim([ymin ymax])

if plot_scales
    if strcmp(test_case,"drake")
        %plot_scale(deltaRh*KM,"\delta_{Rh}");
        %plot_scale(lambda1*KM,"\lambda_1");
        %plot_scale(deltaSM*KM,"\delta_{SM}");
        %plot_scale(deltaM*KM,"\delta_{M}");
    elseif strcmp(test_case,"jet")
        plot_scale(deltaRh*KM,"\delta_{R}");
        plot_scale(lambda1*KM,"\lambda_1");
        plot_scale(deltaSM*KM,"\delta_{M}");
    end
end

set (gca,"fontsize",20);
xlabel("\lambda [km]")

if strcmp(type,"flux")
    ylabel("\Pi(\lambda)/\Pi_{max}")
elseif strcmp(type,"transfer")
    ylabel("T(\lambda)/T_{max}")
else
    ylabel("S(\lambda)")
end
if strcmp(test_case,"drake")
    hp = hp(isgraphics(hp));   % remove invalid handles
    legend(hp, 'Location', 'best', 'FontName', 'Menlo');
end
set (gca,"Xdir","reverse");

if strcmp(test_case,"drake")
    if tanh_strat
        title("Tanh stratification")
    else
        title("Constant/linear stratification")
    end
end

%% Plot power law profile with depth
if drake
    cp1 = 34; cp2 = 34;
    z = -H + H/zlevels * (0.5 + (0:59)); % evenly spaced
else
    cp1 = 271; cp2 = 271;
    z = [-1.25 -3.75 -6.27 -8.82 -11.4 -14.1 -16.8 -19.6 -22.4 -25.4 -28.6 ...
        -31.9 -35.3 -39 -42.8 -47 -51.4 -56.1 -61.3 -66.8 -73 -79.3 -86.4 ...
        -94.2 -103 -112 -122 -134 -146 -160 -176 -193 -212 -233 -257 -283 ...
        -313 -345 -382 -423 -469 -520 -578 -642 -715 -795 -886 -988 -1100 -1230 ...
        -1370 -1540 -1720 -1920 -2150 -2400 -2690 -3010 -3380 -3780];
    z = flipud(z');
end

plot(mean(pow_law(cp1:cp2,:),1),z,'ro-','linewidth',2);grid on;hold on;
axis([-3.5 -1 -H 0])
xlabel('$p, E(k)\propto k^{p}$','interpreter','latex');ylabel('z [m]');set(gca,'fontsize',18)

%% Plot local region
local     = false;
region    = "mid";
data_file = load_data(test_case, dir, run_id, cp_id, zlev, type, machine, false);

fid = fopen(data_file);
data = fread(fid,"double"); 
N = round(-3/2 + sqrt(2*(numel(data)+1)));
data = reshape (data, N+1, N/2+1)';
dmin = min(min(data)); dmax = max(max(data));
fprintf("Minimum value of data = %8.4e\n", dmin);
fprintf("Maximum value of data = %8.4e\n", dmax);

lon  = (-N/2:N/2) * 360/N;
lat  = (-N/4:N/4) * 360/N;

% Target area for local spectral analysis
if local
    if strcmp(region,"southern") % Vortical region at southern edge of land mass
        lat0   = -40;
        lon0   =  20;
        theta0 =  30;
    elseif strcmp(region,"equator") % Vortical region at equator
        lat0   = 0;
        lon0   = 35;
        theta0 = 20;
    elseif strcmp(region,"mid") % Vortical region at 45 N
        lat0   =  45;
        lon0   =  35;
        theta0 =  20;
    elseif strcmp(region,"laminar") % Laminar region
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

field_limit = max(abs(dmin),abs(dmax));

% Color vorticity beyond vort_limit to maximum colors
data(data > field_limit)  =  field_limit;
data(data < -field_limit) = -field_limit;

c_scale = linspace(-field_limit, field_limit, 100); 
smooth  = 0;
lines   = 0;
shift   = 0;
v_title = strip(type,'right','u')+" \bfu";
plot_lon_lat_data(data(1:4:end,1:4:end), lon(1:4:end), lat(1:4:end), c_scale, v_title, smooth, shift, lines);
axis(ax);
%%
function powerlaw (scales, power, range, p)
% Plots power law -p between scales s1 and s2 (s2 > s1) 
[~,k1] = min(abs(scales-range(1))); [~,k2] = min(abs(scales-range(2))); knorm = round((k1+k2)/2);

if abs(p+5/3) < 1e-2
    str = "k^{"+ "-"+ 5 +"/" + 3 +"}";
elseif mod(p,1) ~= 0
    str = "k^{" + compose("%1.2f",p) + "}";
else
    str = "k^{" + compose("%1.0f",p) + "}";
end

loglog(scales(k1:k2),scales(k1:k2).^(-p) * power(knorm)/scales(knorm)^(-p),...
    'k', "linewidth", 2, "DisplayName", str);
end

function plot_scale (scale,name)
y = ylim;
x = scale;
h = loglog([x x], y, "linewidth",1.5,'k');
set(get(get(h,"Annotation"),"LegendInformation"),"IconDisplayStyle","off");
text(0.92*scale, 2.5*y(1), name, "fontsize", 16)
end

function [H, lambda0,lambda1, deltaS, deltaSM, deltaRh, deltaM, radius, Z] = params()
% Physical parameters of simulation

global KM; global tanh_strat

if tanh_strat
    H_linear = 3800;     % depth of linear scaling range (constant/linear)
else
    H_linear = 300;      % depth of linear scaling range (tanh)
end

H_mixed     =  200;     % depth of mixed layer
H_mode      = 1000;     %

uwbc        =  1.0;     % velocity scale

Laplace     =  2;       % 1 = Laplacian, 2 = bi-Laplacian
C_visc      =  1e-3;    % non-dimensional viscosity
dx          =  5e3;     % minimum grid size
dt          =  674;     % time step
g           =  9.80616;
drho        = -4;
ref_density =  1030;
H           =  4e3;     % full depth

visc        =  8.87e8;
scale_omega =  6;
scale_earth =  6;
omega       =  7.29211e-5/scale_omega;
radius      =  6371.229e3/scale_earth;
theta       =  40; % latitude at which to calculate f0 and beta
f0          =  2*omega*sin(deg2rad(theta));
beta        =  2*omega*cos(deg2rad(theta))/radius;
r_b         =  4e-4; % bottom friction

% Brunt-Vaisala frequency
N_bv        = sqrt (- g * drho/ref_density / H_linear);

c0          = sqrt(g*H);
c1          = N_bv * H_mode / pi;
deltaM      = (visc/beta)^(1/(2*Laplace+1)); % Munk layer

% Lengthscales
KM = 1e-3;
lambda0    = c0/f0;             % external radius of deformation
lambda1    = c1/f0;             % internal radius of deformation
deltaS     = r_b/beta;          % Stommel layer
deltaSM    = uwbc/f0;           % submesoscale
deltaRh    = sqrt(uwbc/beta); % Rhine layer
Rey        = uwbc * deltaSM^(2*Laplace-1)/visc; % Reynolds number
Ro         = uwbc / (deltaM*f0); % Rossby number
fprintf('\nlambda0 = %2.1f km lambda1 = %2.1f km\n\n',lambda0*KM,lambda1*KM)
fprintf('deltaS = %2.1f km deltaI = %2.1f km deltaM = %2.1f km deltaSM = %2.1f km\n\n',...
    deltaS*KM, deltaRh*KM, deltaM*KM, deltaSM*KM)
fprintf('Rey = %2.2e Ro = %2.2e N_bv = %2.2e\n\n', Rey, Ro, N_bv)

% Layer depths
Z_tanh = [ ...
    -3.8000e+03  -3.4300e+03  -3.1000e+03  -2.8000e+03  -2.5400e+03  -2.3000e+03  -2.0900e+03  -1.9000e+03  -1.7300e+03  -1.5800e+03  -1.4400e+03  ...
    -1.3200e+03  -1.2100e+03  -1.1100e+03  -1.0200e+03  -9.4600e+02  -8.7400e+02  -8.0900e+02  -7.5100e+02  -6.9800e+02  -6.5000e+02  -6.0600e+02  ...
    -5.6600e+02  -5.3000e+02  -4.9700e+02  -4.6600e+02  -4.3800e+02  -4.1200e+02  -3.8900e+02  -3.6600e+02  -3.4600e+02  -3.2600e+02  -3.0800e+02  ...
    -2.9100e+02  -2.7500e+02  -2.6000e+02  -2.4500e+02  -2.3100e+02  -2.1800e+02  -2.0500e+02  -1.9300e+02  -1.8100e+02  -1.7000e+02  -1.5900e+02  ...
    -1.4800e+02  -1.3800e+02  -1.2700e+02  -1.1700e+02  -1.0700e+02  -9.7400e+01  -8.7800e+01  -7.8300e+01  -6.8800e+01  -5.9500e+01  -5.0200e+01  ...
    -4.1000e+01 -3.1900e+01  -2.2700e+01  -1.3600e+01  -4.5400e+00];


Z_linear = [ ...
    -3.7900e+03  -3.3900e+03  -3.0400e+03  -2.7200e+03  -2.4400e+03  -2.1900e+03  -1.9600e+03  -1.7600e+03  -1.5800e+03  -1.4200e+03  -1.2800e+03  ...
    -1.1600e+03  -1.0400e+03  -9.4300e+02  -8.5300e+02  -7.7200e+02  -7.0000e+02  -6.3600e+02  -5.7800e+02  -5.2600e+02  -4.8000e+02  -4.3800e+02  ...
    -4.0100e+02  -3.6700e+02  -3.3700e+02  -3.0900e+02  -2.8500e+02  -2.6200e+02  -2.4200e+02  -2.2300e+02  -2.0700e+02  -1.9100e+02  -1.7700e+02  ...
    -1.6500e+02  -1.5300e+02  -1.4200e+02  -1.3200e+02  -1.2300e+02  -1.1400e+02  -1.0600e+02  -9.8400e+01  -9.1300e+01  -8.4600e+01  -7.8300e+01  ...
    -7.2300e+01  -6.6500e+01  -6.1000e+01  -5.5800e+01  -5.0700e+01  -4.5800e+01  -4.1000e+01  -3.6400e+01  -3.1900e+01  -2.7500e+01  -2.3100e+01  ...
    -1.8800e+01  -1.4600e+01  -1.0400e+01  -6.2200e+00  -2.0700e+00];

if tanh_strat
    Z = Z_tanh;
else
    Z = Z_linear;
end

end





