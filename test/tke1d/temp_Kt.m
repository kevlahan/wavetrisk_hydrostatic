%% TKE profile
clear all

kato = 0;
time = 1;

remote_file = ['~/hydro/tke/tke1d.6.' sprintf( '%04d', time)];
local_file  = remote_file;
scp_cmd     = ['scp if.mcmaster.ca:' remote_file ' ' local_file];
unix (sprintf(scp_cmd));
T = load(local_file);

% Theoretical value of mixing depth
u_s = 0.01; N_0 = 0.01; t = 30;

if kato
    Dm = - 1.05 * u_s * sqrt(t*60^2/N_0); % Kato and Phillips (1969)
else
    Dm = -11.5;                           % Willis and Deardorff (1974)
end

%figure; 
figure; sgtitle("\Delta t = 1200s N = 100 levels",'fontsize', 24);

% Use Florian's dimensions

subplot(1,2,1);plot(T(:,6),T(:,5),'b-','linewidth',2);hold on;
plot(T(:,6),ones(size(T(:,6)))*Dm,'r--','linewidth',2);
xlabel('Temperature [\circC]');ylabel('Depth [m]'); set(gca,'fontsize',16); grid on;
if kato
    axis([13.25 15.4 -52.4 2.4]); % Kato and Phillips (1969)
else
    axis([19.5 21.0 -25 0]); % Willis and Deardorff (1974)
    %axis([16 20.0 -50 0]); % Willis and Deardorff (1974)
end
pbaspect([0.64 1 1]);

subplot(1,2,2);plot(T(:,3),T(:,2),'b-','linewidth',2);
xlabel('Kt [m^2/s]');ylabel('Depth [m]'); set(gca,'fontsize',16); grid on;
if kato
    axis([-1e-3 2.25e-2 -52.4 2.4]); % Kato and Phillips (1969)
else
    axis([-1e-3 2.5e-2 -25 0]); % Willis and Deardorff (1974)
end
pbaspect([0.64 1 1]);

% subplot(1,3,3);plot(T(:,4),T(:,2),'b-','linewidth',2);
% xlabel('Kv [m^2/s]');ylabel('Depth [m]'); set(gca,'fontsize',16); grid
% axis([0 2e-2 -50 1.2]);
% pbaspect([0.6 1 1]);