%% Load data
clear
%figure
%machine  = "if.mcmaster.ca";
machine   = "nia-datamover1.scinet.utoronto.ca";

test_case="drake"; dir="~/hydro/drake";

% Transfer all spectrum files at once
scp_cmd = "scp "+machine+":"""+dir+"/*spec"" "+dir;
if ~strcmp(machine,"mac")
    unix (sprintf(scp_cmd));
end

%% Analyze spectrum data
type      = "curlu"; 
run_id    = "60layer";
cp_min    = 30;
cp_max    = 30;
zmin      = 59;
zmax      = 59;
plot_spec = false;     % plot spectrum
power     = true;     % plot power law fit
avg       = false;    % plot averaged spectrum
col_spec  = "b-";     % colour for energy spectrum
col_power = "r-";     % colour for power law
%range     = [200 55]; % range for power law fit
%range     = [430 150]; % range for power law fit

% Set physical parameters
[lambda0,lambda1, deltaS, deltaSM, deltaI, radius] = params(test_case)

for cp_id = cp_min:cp_max
    for zlev = zmin:zmax
        name_type = "Layer "+zlev;
        cp        = compose("%04d",cp_id);
        k         = compose("%04d",zlev);

        % Load spectrum data
        if avg % average spectrum
            file_base = run_id+"_"+k+"_"+type;
        else
            file_base = run_id+"_"+cp+"_" +k+"_"+type;
        end
        spec_file  = "~/hydro/"+test_case+"/"+file_base+"_spec";
        pspec = load(spec_file);

        % Plot energy spectra
        if strcmp(type,"u") % velocity spectrum
            pspec(:,2) = pspec(:,2);
        else                % convert vorticity spectrum to energy spectrum integrated over shells
            pspec(:,2) = pspec(:,2)./pspec(:,1).^2;
        end
        scales = 2*pi*radius/1e3./sqrt(pspec(:,1).*(pspec(:,1)+1)); % equivalent length scale (Jeans relation)

        % Fit power law
        if ~exist('range','var') % use default range
            if strcmp(test_case,"drake")
                range = [lambda1 deltaSM]; col_power = "r-"; % colour for power law
            elseif strcmp(test_case,"jet")
                range = [deltaI lambda1]; col_power = "r-"; % colour for power law
            end
        end

        fit_indices = find(scales > range(2) & scales < range(1));
        [P,S] = polyfit(log10(scales(fit_indices)),log10(pspec(fit_indices,2)),1);

        st_err = sqrt(diag(inv(S.R)*inv(S.R'))*S.normr^2/S.df); % error in coefficients from covariance matrix of P

        fprintf("\nFitted power law for checkpoint %d at zlevel %d is %.2f +/- %.2f\n",...
            cp_id, zlev, -P(1), st_err(1));

        pow_law(cp_id,zlev) = -P(1);

        if plot_spec
            loglog(scales,pspec(:,2),col_spec,"linewidth",3,"DisplayName",name_type);hold on;grid on;

            axis([4e0 1e4 1e-10 1e0]);
            set (gca,"fontsize",20);
            xlabel("\lambda (km)");ylabel("S(\lambda)");
            set (gca,"Xdir","reverse");legend;

            if strcmp(test_case,"drake")
                plot_scale(lambda1,"\lambda_1");
                plot_scale(deltaSM,"\delta_{SM}");
            elseif strcmp(test_case,"jet")
                plot_scale(deltaI,"\delta_{I}");
                plot_scale(lambda1,"\lambda_1");
                plot_scale(deltaSM,"\delta_{M}");
            end
            % Plot fit
            if power
                powerlaw (scales, 1.5*pspec(:,2), range, -P(1), col_power)
            end
        end
    end
end
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
function powerlaw (scales, power, range, p, col)
% Plots power law -p between scales s1 and s2 (s2 > s1) with color col
[~,k1] = min(abs(scales-range(1))); [~,k2] = min(abs(scales-range(2))); knorm = round((k1+k2)/2);

if abs(p+5/3) < 1e-2
    str = "k^{"+ "-"+ 5 +"/" + 3 +"}";
elseif mod(p,1) ~= 0
    str = "k^{" + compose("%1.2f",p) + "}";
else
    str = "k^{" + compose("%1.0f",p) + "}";
end

loglog(scales(k1:k2),scales(k1:k2).^(-p) * power(knorm)/scales(knorm)^(-p),col,"linewidth",3,"DisplayName",str);
end

function plot_scale (scale,name)
y = ylim;
x = scale;
h=loglog([x x], y, "k","linewidth",1.5);
set(get(get(h,"Annotation"),"LegendInformation"),"IconDisplayStyle","off");
text(0.92*scale,10*y(1),name,"fontsize",16)
end

function [lambda0,lambda1, deltaS, deltaSM, deltaI, radius] = params(test_case)
% Physical parameters of simulation

if strcmp(test_case,"drake")
    visc        =  99;
    uwbc        =  1.5;
    scale_omega =  6;
    scale_earth =  6;
    omega       =  7.29211e-5/scale_omega;
    radius      =  6371.229e3/scale_earth;
    g           =  9.80616;
    drho        = -4;
    ref_density =  1030;
    H1          =  3e3;
    H2          =  1e3;
    H           =  H1 + H2;
    theta       =  45; % latitude at which to calculate f0 and beta
    f0          =  2*omega*sin(deg2rad(theta));
    beta        =  2*omega*cos(deg2rad(theta))/radius;
    %r_b         =  1.3e-8; % two-layer
    r_b         =  4e-4;
    c0          =  sqrt(g*H);
    %c1          =  sqrt (g*abs(drho)/ref_density * H2*(H-H2)/H); % two-layer
    c1          =  4.8; % m/s
    deltaM      = (visc/beta)^(1/3)/1e3; % Munk layer
elseif strcmp(test_case,"jet")
    visc        =  1.63e7; % hyperviscosity
    uwbc        =  1.4;
    omega       =  1e-4;
    radius      =  1000e3;
    g           =  9.80616;
    drho        = -4;
    ref_density =  1027.8;
    H           =  4e3;
    theta       =  45; % latitude at which to calculate f0 and beta
    f0          =  2*omega*sin(deg2rad(theta));
    beta        =  2*omega*cos(deg2rad(theta))/radius;
    r_b         =  5e-3;
    c0          =  sqrt(g*H);
    c1          =  3.16; % m/s
    deltaM      = (visc/beta)^(1/5)/1e3; % Munk layer
end

% Lengthscales (km)
lambda0    = c0/f0/1e3;             % external radius of deformation
lambda1    = c1/f0/1e3;             % internal radius of deformation
deltaS     = r_b/beta/1e3;        % Stommel layer
deltaSM    = uwbc/f0/1e3;           % submesoscale
deltaI     = sqrt(uwbc/beta)/1e3;   % inertial layer
end



