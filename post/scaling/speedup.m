%
% Parallel scaling results for bbdserv and niagara
%
% Held-Suarez climate test case: J8, 32 layers. No I/O, no topo.
clear; close all; clc
J        = 8;                 % scale
n_col    = 10*4^J + 2;        % number of columns
nz       = 32;                % number of vertical layers 
n_grid   = n_col * nz         % number of grid points

%
% 2560 domains
%

% bbserv
n_dt_bb       = 25; 
cores_bb      = [2 5 8 10 20 40 64 80 100 120 128 140 160 192 200 224 240 256 280 300 320 384]';
cpu_bb        = [2509.5 1267.3 942.4 781.6 456.5 215.6 156.6 96.5 75.3 65.0 ...
    67.7 56.0 53.3 42.5 35.1 34.1 28.0 26.1 23.1 21.0 18.4 15.8]';
% Linearly extrapolate cpu for one core
cpu_bb        = [cpu_bb(1) + (cpu_bb(2)-cpu_bb(1))/(cores_bb(2)-cores_bb(1)) * (1 - cores_bb(1)) ; cpu_bb];
cores_bb      = [1 ; cores_bb];                           
speed_up_bb   = speed_up (cpu_bb);               
abs_perf_bb   = abs_perf (cpu_bb, cores_bb, n_grid, n_dt_bb); 

% niagara
n_dt_nia      = 10;
cores_nia     = [1 2 4 8 10 20 40 80 160 320 640]';
cpu_nia       = [2452.40 1865.0 1036.6 558.4 470.3 249.8 152.6 77.3 43.0 21.5 9.88]';
speed_up_nia  = speed_up (cpu_nia);               
abs_perf_nia  = abs_perf (cpu_nia, cores_nia, n_grid, n_dt_nia); 

% orme
n_dt_orme     = 10;
cores_orme    = [5 10 20 40 60]';
cpu_orme      = [1718.6 937.8 509.2 429.2 263.7]';
speed_up_orme = speed_up (cpu_orme);               
abs_perf_orme = abs_perf (cpu_orme, cores_orme, n_grid, n_dt_orme); 

% 2018 runs (non-adaptive and adaptive)
nz            = 30;         % number of layers
ngrid_nia_old = n_col * nz; % number of grid points
n_dt_nia_old  = 10;         % number of time steps

% non-adaptive
cores_nia_old          = [1 2 4 10 20 40 80 160 320 640 1280 2560]';
cpu_nia_old            = [2.8212E+03 2.3583E+03 1.4431E+03 6.5801E+02 3.8123E+02 1.9202E+02 1.0019E+02 ...
    4.2737E+01 2.1890E+01 1.0355E+01 5.3954E+00 3.0428E+00]';
speed_up_nia_old       = speed_up (cpu_nia_old);               
abs_perf_nia_old       = abs_perf (cpu_nia_old, cores_nia_old, ngrid_nia_old, n_dt_nia_old); 

% adaptive
n_col_nia_adapt_old    = 13762560/4.5/4; % average value
nz_nia_adapt_old       = 18;
ngrid_nia_adapt_old    = n_col_nia_adapt_old * nz_nia_adapt_old;
n_dt_nia_adapt_old     = 300;
cores_nia_adapt_old    = [10 20 40 80 160 320 480 640]';
cpu_nia_adapt_old      = [11806 7458 4995.6 2416.1 1575.4 1016.3 787.01 355.24]';
speed_up_nia_adapt_old = speed_up (cpu_nia_adapt_old);               
abs_perf_nia_adapt_old = abs_perf (cpu_nia_adapt_old, cores_nia_adapt_old, ngrid_nia_adapt_old, n_dt_nia_adapt_old); 

figure
sgtitle({'Strong scaling performance of wavetrisk for Held-Suarez (1994) 0.25^\circ resolution climate',...
    '2560 domains'})
subplot(1,2,1)
loglog(cores_nia,           speed_up_nia,           'bo', 'MarkerSize', 10); hold on
loglog(cores_nia_old,       speed_up_nia_old,       'b*', 'MarkerSize', 10) 
loglog(cores_nia_adapt_old, speed_up_nia_adapt_old, 'md', 'MarkerSize', 10) 
loglog(cores_orme,          speed_up_orme,          'ks', 'MarkerSize', 10) 
loglog(cores_bb,            speed_up_bb,            'ro', 'MarkerSize', 10) 

loglog(cores_nia_old, cores_nia_old/2.5, 'g-')
loglog(cores_nia_old, cores_nia_old/6,   'g-')
loglog(cores_nia_old, cores_nia_old/18,  'g-')

%loglog(cores_nia(5:end),cores_nia(5:end)/2.6,'g-')
legend('niagara','2018 niagara','2018 niagara adapt', 'orme', 'bbserv','perfect')
xlabel('Number of cores'); ylabel('Speed up'); grid on; 
axis([1e0 1e4 1e0 1e4])
set(gca,'fontsize',16)

subplot(1,2,2)
loglog(cores_nia,           abs_perf_nia,           'bo-', 'MarkerSize', 10); hold on
loglog(cores_nia_old,       abs_perf_nia_old,       'b*-', 'MarkerSize', 10) 
loglog(cores_nia_adapt_old, abs_perf_nia_adapt_old, 'md-', 'MarkerSize', 10)
loglog(cores_orme,          abs_perf_orme,          'ks-', 'MarkerSize', 10)
loglog(cores_bb,            abs_perf_bb,            'ro-', 'MarkerSize', 10)
legend('niagara','2018 niagara','2018 niagara adapt', 'orme', 'bbserv')
xlabel('Number of cores'); ylabel('Absolute performance (ms/dt/cell)'); grid on;
axis([1e0 1e4 1e-3 1e-1])
set(gca,'fontsize',16)

%
% 640 domains 
%

% bbserv
n_dt_bb      = 25;   
cores_bb     = [2 5 8 10 20 40 64 80 100 120 128 140 160 192 200 224 240 256 280 300 320 384]';
cpu_bb       = [1032.4 440.6 291.9 233.4 119.4 57.6 41.4 30.4 27.6 22.4 19.4 19.0 15.7 15.3 15.2 ...
     12.5 12.3 12.4 11.9 11.8 9.12 8.31]';
% Linearly extrapolate cpu for one core
cpu_bb       = [cpu_bb(1) + (cpu_bb(2)-cpu_bb(1))/(cores_bb(2)-cores_bb(1)) * (1 - cores_bb(1)) ; cpu_bb];
cores_bb     = [1 ; cores_bb];
speed_up_bb  = speed_up (cpu_bb);               
abs_perf_bb  = abs_perf (cpu_bb, cores_bb, n_grid, n_dt_bb); 

% niagara
n_dt_nia      = 10;
cores_nia     = [1 2 4 8 10 20 40 80 160 320 640]';
cpu_nia       = [1136.5 615.6 307.6 158.8 128.3 65.7 43.8 24.0 12.9 7.49 4.06]';
speed_up_nia  = speed_up (cpu_nia);               
abs_perf_nia  = abs_perf (cpu_nia, cores_nia, n_grid, n_dt_nia); 

% orme
n_dt_orme     = 10;
cores_orme    = [4 10 20 40 60]';
cpu_orme      = [1657.8 847.8 447.6 263.8 228.6]';speed_up_orme = speed_up (cpu_orme);               
abs_perf_orme = abs_perf (cpu_orme, cores_orme, n_grid, n_dt_orme); 

figure
sgtitle({'Strong scaling performance of wavetrisk for Held-Suarez (1994) 0.25^\circ resolution climate',...
    '640 domains'})
subplot(1,2,1)
loglog(cores_nia, speed_up_nia,    'bo', 'MarkerSize', 10); hold on
loglog(cores_orme,  speed_up_orme, 'ks', 'MarkerSize', 10); 
loglog(cores_bb,  speed_up_bb,     'ro', 'MarkerSize', 10);

loglog(cores_nia, cores_nia/1.1, 'g-')
loglog(cores_bb,  cores_bb/5,    'g-')
loglog(cores_bb,  cores_bb/1.9,  'g-')

legend('niagara','orme','bbserv','perfect')
xlabel('Number of cores'); ylabel('Speed up'); grid on;
axis([1e0 1e3 1e0 1e3])
set(gca,'fontsize',16)

subplot(1,2,2)
loglog(cores_nia,  abs_perf_nia,  'bo', 'MarkerSize', 10); hold on
loglog(cores_orme, abs_perf_orme, 'ks', 'MarkerSize', 10);
loglog(cores_bb,   abs_perf_bb,   'ro', 'MarkerSize', 10);
legend('niagara','orme','bbserv')
xlabel('Number of cores'); ylabel('Absolute performance [ms]'); grid on; 
axis([1e0 1e3 1e-3 1e-1])
set(gca,'fontsize',16)

%% Results from other codes (varies from 5e-5 for 8 cores to 1.6e-1 for 32768 cores)
clear; clc
day  = 24*60^2;    % day in seconds

% Lupo etal (2017)
clear; clc
day   = 24*60^2;    % day in seconds
dt    = 30;         % time step in seconds
n_dt  = 5*day/dt;   % number of time steps in the simulation
n_grid = 321*640*32; % number of grid points

fprintf('\n%s\n\n','Absolute performance from other codes [ms/dt/cell]')
fprintf('%s %.2e %s\n','Lupo et al (2017):              ', abs_perf (42174, 256, n_grid, n_dt), '(256 cores)')
fprintf('%s %.2e %s\n','                                ', abs_perf (30337, 128, n_grid, n_dt), '(128 cores)')
fprintf('%s %.2e %s\n','                                ', abs_perf (26695,  64, n_grid, n_dt), '( 64 cores)')
fprintf('%s %.2e %s\n','                                ', abs_perf (22029,  32, n_grid, n_dt), '( 32 cores)')
fprintf('%s %.2e %s',  '                                ', abs_perf (28959,  16, n_grid, n_dt), '( 16 cores)')

% Zuo et al (2007) benchmark 3
n_dt  = 200;
n_grid = 2048*256*30;
cpu   = 20;
cores = 8;
fprintf ('\n\n%s %.2e %s','Zuo et al (2007) benchmark 3:   ',abs_perf(cpu, cores, n_grid, n_dt), '(8 cores)')

% Wang et al (2005)
n_grid     = 1024^2*20;
cpu   = 327;
cores = 64;
n_dt  = 6000;
fprintf ('\n\n%s %.2e %s\n','Wong et al (2005):              ',abs_perf(cpu, cores, n_grid, n_dt), '(64 cores)')

% Heinzeller
NZ = 41;

% 2 km
cpu = 3585;
n_grid = 147456002*NZ;
cores = 2048*16;
dt = 5;
n_dt = 10*60/dt;
fprintf ('\n%s %.2e %s','Heinzeller:                     ',abs_perf(cpu,  cores, n_grid, n_dt), '(32768 cores,  2 km)')

% 120 km
n_grid = 40962*NZ;
dt = 150;
n_dt = day/dt;
cores = 160;
cpu   = 236;
fprintf ('\n%s %.2e %s\n','                                ',abs_perf(cpu,  cores, n_grid, n_dt), '(160 cores, 120 km)')

% Jablonowski and Williamson (2006)
core = 32;         % number of cores
NZ   = 26;         % number of vertical levels

% CAM3 Euler
cpu = 483;
dt  = 300; 
n_grid   = 256*512*NZ;
fprintf ('\n%s %.2e %s\n','Jablonowski & Williamson (2006):', abs_perf(cpu,  cores, n_grid, day/dt), '(32 cores, CAM3 Euler)')

% CAM3 finite volume
cpu = 625;
dt  = 90;
n_grid   = 361*576*NZ;
fprintf ('%s %.2e %s\n','                                ', abs_perf(cpu,  cores, n_grid, day/dt), '(32 cores, CAM3 FV)')

% GME 
cpu = 325;
dt  = 200;
n_grid   = 163842*NZ; 
fprintf ('%s %.2e %s\n\n','                                ', abs_perf(cpu, cores, n_grid, day/dt), '(32 cores, GME)')


%%% Functions

function perf = abs_perf (cpu, cores, n_grid, n_dt)
% Absolute performance (ms/dt)
perf = (cpu/n_dt) ./ (n_grid./cores) * 1e3;
end

function speed = speed_up (cpu)
% Speed up
speed = cpu(1)./cpu;
end
