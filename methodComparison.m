% load lattice data
load('results_2d_square_2017_05_20_12_21_17.mat')
m = [ghParams.m];
lattice_pathLen = 1./m(2:2:end);
lattice_Deff = [results_homog.Deff];
lattice_Deff = lattice_Deff(2:2:end)./lattice_Deff(1:2:end);

% load lattice diag jump data
load('results_2d_square_diagJumps_2017_05_20_12_22_53.mat')
m = [ghParams.m];
latticeD_pathLen = .5*(1./m(2:2:end) + sqrt(2)./m(2:2:end));
latticeD_Deff = [results_homog.Deff];
latticeD_Deff = latticeD_Deff(2:2:end)./latticeD_Deff(1:2:end);

% continuous data
load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
specs = [results.sim_spec];
N = [specs.sims];
cont_pathLen = [specs.path_len];
cont_Deff = [results.Deff];
err = 2*sqrt([results.DeffVar]./N);

figure
hold on
plot(log(1./lattice_pathLen),lattice_Deff,'.-','markersize',15)
plot(log(1./latticeD_pathLen),latticeD_Deff,'.-','markersize',15)
%plot(log(1./cont_pathLen),cont_Deff,'.-','markersize',15)
errorbar(log(1./cont_pathLen),cont_Deff,err,'.-','markersize',15);

xlabel('log(1/path length)')
ylabel('Deff')
legend('Lattice','Lattice (diag jumps)','Continuous (2 std dev)','location','southwest')