% Scaling results: J8, 32 layers, 2560 domains. No I/O, no topo.
clear; close all;
J      = 8;                 % scale
n_col  = 10*4^J + 2;        % number of columns
nz     = 32;                % number of vertical layers

% bbserv
n_dt_bb = 25; % number of time steps

cores_bb = [2 5 8 10 20 40 64 80 100 120 128 140 160 192 200 224 240 256 280 300 320 384]';
cpu_bb   = [2509.5 1267.3 942.4 781.6 456.5 215.6 156.6 96.5 75.3 65.0 ...
    67.7 56.0 53.3 42.5 35.1 34.1 28.0 26.1 23.1 21.0 18.4 15.8]';

% Linearly extrapolate cpu for one core
cpu_bb = [cpu_bb(1) + (cpu_bb(2)-cpu_bb(1))/(cores_bb(2)-cores_bb(1)) * (1 - cores_bb(1)) ; cpu_bb];

cores_bb = [1 ; cores_bb];

N_bb = nz * n_col./cores_bb; % normalization for absolute performance

speed_up_bb = cpu_bb(1)./cpu_bb;      % speed up
abs_perf_bb = cpu_bb/n_dt_bb ./N_bb * 1e3; % absolute performace in ms  

% niagara
n_dt_nia  = 10; % number of time steps

cores_nia = [1 2 4 8 10 20 40 80 160 320 640]';
cpu_nia   = [2452.40 1865.0 1036.6 558.4 470.3 249.8 152.6 77.3 43.0 21.5 9.88]';

N_nia = nz * n_col./cores_nia; % absolute performance measure

speed_up_nia = cpu_nia(1)./cpu_nia;
abs_perf_nia = (cpu_nia/n_dt_nia) ./ N_nia * 1e3;

figure
sgtitle({'Strong scaling performance of wavetrisk for Held-Suarez (1994) 0.25^\circ resolution climate',...
    '2560 domains'})
subplot(1,2,1)
loglog(cores_nia, speed_up_nia,'bo','MarkerSize',10); hold on
loglog(cores_bb,  speed_up_bb, 'ro','MarkerSize',10); hold on
loglog(cores_nia,cores_nia/2.5,'g-')
%loglog(cores_nia(5:end),cores_nia(5:end)/2.6,'g-')
legend('niagara','bbserv','perfect')
xlabel('Number of cores'); ylabel('Speed up'); grid on;
axis([1e0 1e3 1e0 1e3])
set(gca,'fontsize',16)

subplot(1,2,2)
semilogx(cores_nia, abs_perf_nia,'bo','MarkerSize',10); hold on
semilogx(cores_bb,  abs_perf_bb, 'ro','MarkerSize',10);
legend('niagara','bbserv')
xlabel('Number of cores'); ylabel('Absolute performance (ms)'); grid on;
axis([1e0 1e3 0 0.04])
set(gca,'fontsize',16)

%% 640 domains (bbserv only)
J      = 8;                 % scale
n_col  = 10*4^J + 2;        % number of columns
nz     = 32;                % number of vertical layers

n_dt_bb = 25; % number of time steps
cores_bb = [2 5 8 10 20 40 64 80 100 120 128 140 160 192 200 224 240 256 280 300 320 384]';

cpu_bb = [1032.4 440.6 291.9 233.4 119.4 57.6 41.4 30.4 27.6 22.4 19.4 19.0 15.7 15.3 15.2 ...
     12.5 12.3 12.4 11.9 11.8 9.12 8.31]';

% Linearly extrapolate cpu for one core
cpu_bb = [cpu_bb(1) + (cpu_bb(2)-cpu_bb(1))/(cores_bb(2)-cores_bb(1)) * (1 - cores_bb(1)) ; cpu_bb];
cores_bb = [1 ; cores_bb];

N_bb = nz * n_col./cores_bb; % normalization for absolute performance

speed_up_bb = cpu_bb(1)./cpu_bb;      % speed up
abs_perf_bb = cpu_bb/n_dt_bb ./N_bb * 1e3; % absolute performace in ms  

figure
sgtitle({'Strong scaling performance of wavetrisk for Held-Suarez (1994) 0.25^\circ resolution climate',...
    '640 domains'})
subplot(1,2,1)
loglog(cores_bb,  speed_up_bb, 'ro','MarkerSize',10); hold on
loglog(cores_bb,cores_bb/1.9,'g-')
legend('bbserv','perfect')
xlabel('Number of cores'); ylabel('Speed up'); grid on;
axis([1e0 1e3 1e0 1e3])
set(gca,'fontsize',16)

subplot(1,2,2)
semilogx(cores_bb,  abs_perf_bb, 'ro','MarkerSize',10);
xlabel('Number of cores'); ylabel('Absolute performance (ms)'); grid on;
axis([1e0 1e3 0 0.01])
set(gca,'fontsize',16)


