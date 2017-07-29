
load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
results_cont = results;

load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/results/results_2d_square_2017_07_19_09_37_30.mat');
results_discrete = results_homog;
ghParams_discrete = ghParams;

load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/Results_2d_square_diagJumpsCorrected/results_2017_07_29_13_13_26.mat');
results_discrete_diagC = results_homog;

load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/Results_2d_square_diagJumps/results_2017_07_29_13_28_21.mat');
results_discrete_diag = results_homog;

sim_spec = [results_cont.sim_spec];

pathLen_cont = [sim_spec.path_len];
pathLen_disc = 1./[ghParams_discrete.m];
pathLen_disc_diagC = zeros(1,length(results_discrete_diagC));
pathLen_disc_diag = zeros(1,length(results_discrete_diag));
for i = 1:length(results_discrete_diagC)
    curRes = results_discrete_diagC(i);
    
    pathLen = sqrt(sum(curRes.edgeJumps.^2,2));
    meanPathLen = sum(pathLen.*[curRes.edgeRates])./sum([curRes.edgeRates]);
    pathLen_disc_diagC(i) = meanPathLen;
    
    curRes = results_discrete_diag(i);
    
    pathLen = sqrt(sum(curRes.edgeJumps.^2,2));
    meanPathLen = sum(pathLen.*[curRes.edgeRates])./sum([curRes.edgeRates]);
    pathLen_disc_diag(i) = meanPathLen;
end

CI = [results_cont.Deff95CI];
CIlow = CI(1:2:end);
CIhigh = CI(2:2:end);
Deff_cont = [results_cont.Deff];
Deff_disc = [results_discrete.Deff];
Deff_disc_diagC = [results_discrete_diagC.Deff];
Deff_disc_diag = [results_discrete_diag.Deff];

figure

subplot(1,3,1)
hold on
errorbar(pathLen_cont,Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(pathLen_disc,Deff_disc,'.-','markersize',15);
plot(pathLen_disc_diag,Deff_disc_diag,'.-','markersize',15);
plot(pathLen_disc_diagC,Deff_disc_diagC,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','location','southeast');
xlabel('mean path length')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;

subplot(1,3,2)
hold on
errorbar(1./pathLen_cont,Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(1./pathLen_disc,Deff_disc,'.-','markersize',15);
plot(1./pathLen_disc_diag,Deff_disc_diag,'.-','markersize',15);
plot(1./pathLen_disc_diagC,Deff_disc_diagC,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','location','southeast');
xlabel('1 / mean path length')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;

subplot(1,3,3)
hold on
errorbar(log2(1./pathLen_cont),Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(log2(1./pathLen_disc),Deff_disc,'.-','markersize',15);
plot(log2(1./pathLen_disc_diag),Deff_disc_diag,'.-','markersize',15);
plot(log2(1./pathLen_disc_diagC),Deff_disc_diagC,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','location','southeast');
xlabel('log(1 / mean path length)')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;
ca.XLim(1) = 1;
%{
yca = get(gca,'ylabel');
yca.Position = [ -30.912402524544123 0.720000046491606  -1 ];

xca = get(gca,'xlabel');
xca.Position = [ 149.9891430511475 0.6631063628431 -0.01 ];
%}
%lh.Position = 