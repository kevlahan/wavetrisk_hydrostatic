Table = readtable('SimpleJ5J7Z30_profile.csv');
VarNames= Table.Properties.VariableNames

OMEGA_Idx = find(strcmp(VarNames, 'P_Ps'));
Ps_Idx    = find(strcmp(VarNames, 'Ps'));
P_Ps_Idx  = find(strcmp(VarNames, 'P_Ps'));
T_Idx     = find(strcmp(VarNames, 'Temperature'));
u_Idx     = find(strcmp(VarNames, 'VelocityZonal'));
v_Idx     = find(strcmp(VarNames, 'VelocityMeridional'));
vort_Idx  = find(strcmp(VarNames, 'Vorticity'));
dz_Idx    = find(strcmp(VarNames, 'dz'));

semilogy(A{:,T_Idx},A{:,P_Ps_Idx}*1e5,'r-','LineWidth',1.5)
grid on
set(gca, 'YDir','reverse')
xlabel('Temperature [K]'); ylabel('P/P_S')
