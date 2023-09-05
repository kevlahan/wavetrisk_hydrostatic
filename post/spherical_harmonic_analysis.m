%% Load data
clear; clc
%figure
%machine  = "if.mcmaster.ca";
machine   = "nia-datamover1.scinet.utoronto.ca";

dir_remote = "~/hydro/drake";                dir_local = "~/hydro/drake";
%dir_remote = "~/proj/jet/gmd_paper/spectra"; dir_local = "~/hydro/jet";
%dir_remote = "~/proj/drake/old/J6/2layer_mix_slow"; dir_local = "~/hydro/drake";

% Transfer all spectrum files at once
scp_cmd = "scp "+machine+":"""+dir_remote+"/*spec"" "+dir_local;
if ~strcmp(machine,"mac")
    unix (sprintf(scp_cmd));
end

%% Analyze spectrum data
clear; clc; close all
drake = true;
figure;
if drake
    test_case = "drake";
    run_id    = "12layer";
    cp_min    =  41; cp_max = 41;
    type      = "curlu";
else
    type      = "curlu";
    test_case = "jet";
    run_id    = "jet";
    cp_min    = 271; cp_max = 271;
end

% Set physical parameters
KM = 1e-3;
[H, lambda0,lambda1, deltaS, deltaSM, deltaI, deltaM, radius] = params(test_case);

zlevels   = 12;
zmin      = 7;
zmax      = 7;

plot_spec = true;     % plot spectrum
power     = true;     % plot power law fit
avg       = false;    % plot averaged spectrum
col_spec  = "b-";     % colour for energy spectrum
col_power = "r-";     % colour for power law

range     = [deltaI*KM lambda1*KM]; % range for power law fit
%range     = [lambda1*KM deltaSM*KM]; % range for power law fit

%range     = [28 8]; % range for power law fit

pow_law = zeros(cp_max-cp_min+1,zlevels);
for cp_id = cp_min:cp_max
    for zlev = zmin:zmax
        % Load spectrum data
        name_type = "Layer "+zlev;
        cp        = compose("%04d",cp_id);
        k         = compose("%04d",zlev);
        if avg % average spectrum
            file_base = run_id+"_"+k+"_"+type;
        else
            file_base = run_id+"_"+cp+"_"+k+"_"+type;
            %file_base = run_id+"_"+cp+"_"+type;
        end
        
        spec_file = "~/hydro/"+test_case+"/"+file_base+"_spec";
        try
            pspec = load(spec_file);
        catch ME
            fprintf('\n File %s not present ... continuing\n',spec_file)
            pause
            continue
        end

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
                range = [deltaI lambda1]; col_power = "r-"; % colour for power law
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

            axis([1e0 1e4 1e-10 1e1]);
            set (gca,"fontsize",20);
            xlabel("\lambda (km)");ylabel("S(\lambda)");
            set (gca,"Xdir","reverse");legend;

            if strcmp(test_case,"drake")
                plot_scale(lambda1*KM,"\lambda_1");
                plot_scale(deltaSM*KM,"\delta_{SM}");
                plot_scale(deltaI*KM,"\delta_{I}");
                plot_scale(deltaM*KM,"\delta_{M}");
            elseif strcmp(test_case,"jet")
                plot_scale(deltaI*KM,"\delta_{I}");
                plot_scale(lambda1*KM,"\lambda_1");
                plot_scale(deltaSM*KM,"\delta_{M}");
            end
            % Plot fit
            if power
                powerlaw (scales, 1.5*pspec(:,2), range, -P(1), col_power)
            end
        end
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

function [H, lambda0,lambda1, deltaS, deltaSM, deltaI, deltaM, radius] = params(test_case)
% Physical parameters of simulation

if strcmp(test_case,"drake")
    Laplace     =  2;    % 1 = Laplacian, 2 = bi-Laplacian
    C_visc      =  1e-3; % non-dimensional viscosity
    dx          =  1.25e3;  % minimum grid size
    dt          =  126;  % time step
    visc        =  C_visc * dx^(2*Laplace)/dt;
    uwbc        =  0.8;  % velocity scale
    scale_omega =  1;
    scale_earth =  6;
    omega       =  7.29211e-5/scale_omega;
    radius      =  6371.229e3/scale_earth;
    g           =  9.80616;
    drho        = -2;
    ref_density =  1030;
    H1          =  3e3;
    H2          =  1e3;
    H           =  H1 + H2;
    theta       =  40; % latitude at which to calculate f0 and beta
    f0          =  2*omega*sin(deg2rad(theta));
    beta        =  2*omega*cos(deg2rad(theta))/radius;
    r_b         =  7e-8; % bottom friction
    c0          =  sqrt(g*H);
    c1          =  2.85; % internal wave speed
    deltaM      = (visc/beta)^(1/(2*Laplace+1)); % Munk layer
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

% Lengthscales
KM = 1e-3;
lambda0    = c0/f0;             % external radius of deformation
lambda1    = c1/f0;             % internal radius of deformation
deltaS     = r_b/beta;          % Stommel layer
deltaSM    = uwbc/f0;           % submesoscale
deltaI     = sqrt(uwbc/beta);   % inertial layer
Rey    = uwbc*deltaSM^(2*Laplace-1)/visc;

fprintf('\nlambda0 = %2.1f km lambda1 = %2.1f km\n',lambda0*KM,lambda1*KM)
fprintf('\ndeltaS = %2.1f km deltaI = %2.1f km deltaM = %2.1f km deltaSM = %2.1f km\n',...
    deltaS*KM, deltaI*KM, deltaM*KM, deltaSM*KM)
fprintf('\nRey = %2.2e\n', Rey)
visc^(1/(2*Laplace+1)) * KM
end



