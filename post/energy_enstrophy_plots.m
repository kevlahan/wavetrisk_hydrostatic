clear all;
run_id    = '1layer_J6';
% run_id  = '2layer_J6';
test_case = 'drake';
itype     = 'barotropic vorticity' % field to analyze
machine   = 'niagara.computecanada.ca'

file_base   = [run_id];
remote_file = ['~/hydro/' test_case '/' file_base  '_*en*'];
scp_cmd     = ['scp "' machine ':' remote_file '"  ~/hydro/' test_case '/.'];
unix (sprintf(scp_cmd));
%% 2 layer case
% Plot kinetic energy
e = load(['~/hydro/' test_case '/' file_base  '_kinetic_energy']);
figure(1);plot(e(:,1),e(:,4),'b',e(:,1),100*e(:,5),'r',e(:,1),e(:,2),'k',e(:,1),e(:,3),'g','linewidth',1.4);
legend('Barotropic','100 * Baroclinic','Layer 1','Layer 2');
title('Kinetic energies in two layer shallow water system');
xlabel('Time (days)');ylabel('Kinetic energy density');grid on;
%axis([0 1300 0 700]);
set(gca,'FontSize',16);
%%
savefig('~/hydro/drake/kinetic_energy');print -dpng ~/hydro/drake/energy_2layer;
%%
% Plot potential enstrophy
w = load('~/hydro/drake/2layer_fill_pot_enstrophy');
figure(2);plot(w(:,1),w(:,4),'b',w(:,1),w(:,5),'g',w(:,1),w(:,6),'r','linewidth',1.4);
legend('Barotropic','Baroclinic (layer 1)','Baroclinic (layer 2)');
title('Potential enstrophies in two layer shallow water system');
xlabel('Time (days)');ylabel('Potential enstrophy');grid on; 
set(gca,'FontSize',16);
%%
savefig('~/hydro/drake/enstrophy');print -dpng ~/hydro/drake/enstrophy_2layer;
%% 1 layer case
% Plot kinetic energy
e = load(['~/hydro/' test_case '/' file_base  '_kinetic_energy']);
figure(1);plot(e(:,1),e(:,2),'b','linewidth',1.4);
title('Kinetic energy in one layer shallow water system');
xlabel('Time (days)');ylabel('Kinetic energy density');grid on;
%axis([0 1300 0 700]);
set(gca,'FontSize',16);
%%
savefig('~/hydro/drake/kinetic_energy');print -dpng ~/hydro/drake/energy_1layer;
%%
% Plot potential enstrophy
w = load(['~/hydro/' test_case '/' file_base  '_pot_enstrophy']);
figure(2);plot(w(:,1),w(:,2),'b','linewidth',1.4);
title('Potential enstrophy in one layer shallow water system');
xlabel('Time (days)');ylabel('Potential enstrophy');grid on; 
set(gca,'FontSize',16);
%%
savefig('~/hydro/drake/enstrophy');print -dpng ~/hydro/drake/enstrophy_1layer;