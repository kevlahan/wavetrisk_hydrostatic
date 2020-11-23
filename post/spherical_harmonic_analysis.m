% Analyse spherical harmonic data
clear;
test_case  = 'drake';
run_id     = '1layer_J6';
%run_id     = '2layer_fill';
%type       = 'baroclinic_2';
type       = 'barotropic';
cp_id      = '0027';
local      = false;
machine    = 'if.mcmaster.ca';
%machine    = 'mac';
new_transfer = true;
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
loglog(pspec(:,1),pspec(:,2),'r-','linewidth',1.5,'DisplayName','global');hold on;grid on;

if strcmp(type,'barotropic')
    if strcmp(run_id,'1layer_J8')
        p = -2;
        nspec=numel(pspec(:,1));
        k1 = 80; k2 = round(0.8*nspec); knorm = 200;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'m--','linewidth',1.5,'DisplayName','k^{-2}');
    elseif strcmp(run_id,'1layer_J6')
        p = -3; k1 = 3; k2 = 400; knorm = 80;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'b--','linewidth',1.5,'DisplayName','k^{-3}');
        p = -5; k1 = 400; k2 = 2000; knorm = 1000;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-5}');
    elseif strcmp(run_id,'2layer_J8')
        p = -2.2;
        nspec=numel(pspec(:,1));
        k1 = 50; k2 = round(0.8*nspec); knorm = 200;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'m--','linewidth',1.5,'DisplayName','k^{-2.2}');
    elseif strcmp(run_id,'2layer_J6')
        p = -3.3; k1 = 3; k2 = 500; knorm = 80;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-3.3}');
        p = -6; k1 = 300; k2 = 2500; knorm = 800;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-6}');
    elseif strcmp(run_id,'2layer_fill')
        p = -3.3;
        nspec=numel(pspec(:,1));
        k1 = 3; k2 = 200; knorm = 50;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-3.3}');
        
        p = -6;
        nspec=numel(pspec(:,1));
        k1 = 100; k2 = 2400; knorm = 500;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-6}');
    end
end

if strcmp(type,'baroclinic_1') || strcmp(type,'baroclinic_2')
    if strcmp(run_id,'2layer_J8')
        p = -5/3;
        nspec=numel(pspec(:,1));
        k1 = 40; k2 = round(nspec); knorm = 200;
        loglog(pspec(k1:k2,1),0.9*pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'m--','linewidth',1.5,'DisplayName','k^{-5/3}');
    elseif strcmp(run_id,'2layer_fill')
        p = -5/3;
        nspec=numel(pspec(:,1));
        k1 = 4; k2 = 150; knorm = 50;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'b--','linewidth',1.5,'DisplayName','k^{-5/3}');
        
        p = -6;
        nspec=numel(pspec(:,1));
        k1 = 150; k2 = 2400; knorm = 500;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-6}');
    elseif strcmp(run_id,'2layer_J6')
        p = -5/3; k1 = 3; k2 = 500; knorm = 80;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'b--','linewidth',1.5,'DisplayName','k^{-5/3}');
        p = -5; k1 = 300; k2 = 2500; knorm = 800;
        loglog(pspec(k1:k2,1),pspec(k1:k2,1).^p * pspec(knorm,2)/pspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-5}');
    end
end

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
    
    if strcmp(type,'barotropic')
        if strcmp(run_id,'1layer_J8')
            p = -2;
            k1 = 100; k2 = round(nspec); knorm = 200;
            loglog(lspec(k1:k2,1),lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'m--','linewidth',1.5,'DisplayName','k^{-2}');
        elseif strcmp(run_id,'2layer_fill')
            p = -3.3;
            k1 = 10; k2 = round(0.3*nspec); knorm = 50;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-3.3}');
            p = -6;
            k1 = 200; k2 = 700; knorm = 500;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'k--','linewidth',1.5,'DisplayName','k^{-6}');
        elseif strcmp(run_id,'2layer_J6')
            p = -3.3;
            k1 = 10; k2 = round(0.3*nspec); knorm = 50;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-3.3}');
            p = -6;
            k1 = 200; k2 = 700; knorm = 500;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'k--','linewidth',1.5,'DisplayName','k^{-6}');
        end
    end
    
    if strcmp(type,'baroclinic_1')||strcmp(type,'baroclinic_2')
        if strcmp(run_id,'1layer_J8')
            p = -3.3;
            k1 = 40; k2 = round(nspec);knorm=50;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-3.3}');
        elseif strcmp(run_id,'2layer_fill')
            p = -1.2;
            k1 = 3; k2 = 150; knorm = 50;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-1.2}');
            p = -6;
            k1 = 200; k2 = 2000; knorm = 500;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-6}');
        elseif strcmp(run_id,'2layer_J6')
            p = -1.2;
            k1 = 3; k2 = 150; knorm = 50;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'r--','linewidth',1.5,'DisplayName','k^{-1.2}');
            p = -6;
            k1 = 200; k2 = 2000; knorm = 500;
            loglog(lspec(k1:k2,1),1.0*lspec(k1:k2,1).^p * lspec(knorm,2)/lspec(knorm,1)^p,'g--','linewidth',1.5,'DisplayName','k^{-6}');
        end
    end
end

titl = [type ' spherical harmonics spectrum'];
xlabel('l');ylabel('power');title(titl);legend;
axis([1 1e4 1e-12 1e1])
set (gca,'fontsize',18);
