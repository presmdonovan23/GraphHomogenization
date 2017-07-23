load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
results_cont = results;

load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/results/results_2d_square_2017_07_19_09_37_30.mat');
results_discrete = results_homog;

sim_spec = [results.sim_spec];

figure
hold on
plot(1./[sim_spec.path_len],[results_cont.Deff],'.-','markersize',10);
plot([ghParams.m],[results_discrete.Deff],'.-','markersize',10);

lh = legend('Continuous','Discrete','location','southeast');
xlabel('m')
ylabel('Deff')

ca = gca;
ca.FontSize = 10;

yca = get(gca,'ylabel');
yca.Position = [ -30.912402524544123 0.720000046491606  -1 ];

xca = get(gca,'xlabel');
xca.Position = [ 149.9891430511475 0.6631063628431 -0.01 ];

%lh.Position = 