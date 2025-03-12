% Plots vertical profiles of averaged variables produced by lonlat_to_3D.py

Data = readtable('SimpleJ5J7Z30_profile.csv');
VarNames= Data.Properties.VariableNames

OMEGA_Idx = find(strcmp(VarNames, 'P_Ps'));
Ps_Idx    = find(strcmp(VarNames, 'Ps'));
P_Ps_Idx  = find(strcmp(VarNames, 'P_Ps'));
T_Idx     = find(strcmp(VarNames, 'Temperature'));
u_Idx     = find(strcmp(VarNames, 'VelocityZonal'));
v_Idx     = find(strcmp(VarNames, 'VelocityMeridional'));
vort_Idx  = find(strcmp(VarNames, 'Vorticity'));
dz_Idx    = find(strcmp(VarNames, 'dz'));

semilogy(Data{:,T_Idx},Data{:,P_Ps_Idx}.*Data{:,Ps_Idx},'r-','LineWidth',2)
grid on
set(gca, 'YDir','reverse')
xlabel('Temperature [K]'); ylabel('P [Pa]')
axis([120 340 1e3 1e5])
set(gca,'fontsize',18)