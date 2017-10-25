%% Continuous vs discrete vs discrete w/ diag vs discrete w/ diag corrected
m = 2.^[2:8];
m_ext = 2.^[2:9];

load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
results_cont = results;

load('Results/results_2d_square_2017_07_19_09_37_30.mat','results_homog');
results_discrete = results_homog;

load('Results_2d_square_diagJumpsCorrected/results_2017_08_13_18_38_51.mat','results_homog');
results_discrete_diagC = results_homog;

load('Results_2d_square_diagJumps/results_2017_07_29_13_28_21.mat','results_homog');
results_discrete_diag = results_homog;

load('Results_2d_squareBdySlow/results_2017_08_14_12_49_25.mat','results_homog');
results_obsBdy = results_homog;

load('Results_2d_squareBdySlow_diagJumps/results_2017_08_14_12_49_38.mat','results_homog');
results_obsBdy_diag = results_homog;

sim_spec = [results_cont.sim_spec];

pathLen_cont = [sim_spec.path_len];
pathLen_disc = 1./m;
pathLen_disc_diagC = (.5*(sqrt(2) + 1))./m_ext;
pathLen_disc_diag = (.5*(sqrt(2) + 1))./m;
pathLen_obsBdy = 1./m_ext;
pathLen_obsBdy_diag = (.5*(sqrt(2) + 1))./m_ext;

CI = [results_cont.Deff95CI];
CIlow = CI(1:2:end);
CIhigh = CI(2:2:end);
Deff_cont = [results_cont.Deff];
Deff_disc = [results_discrete.Deff];
Deff_disc_diagC = [results_discrete_diagC.Deff];
Deff_disc_diag = [results_discrete_diag.Deff];
Deff_obsBdy = [results_obsBdy.Deff];
Deff_obsBdy_diag = [results_obsBdy_diag.Deff];

figure

hold on
errorbar(log2(1./pathLen_cont),Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(log2(1./pathLen_disc),Deff_disc,'.-','markersize',15);
plot(log2(1./pathLen_disc_diag),Deff_disc_diag,'.-','markersize',15);
plot(log2(1./pathLen_disc_diagC),Deff_disc_diagC,'.-','markersize',15);
plot(log2(1./pathLen_obsBdy),Deff_obsBdy,'.-','markersize',15);
plot(log2(1./pathLen_obsBdy_diag),Deff_obsBdy_diag,'.-','markersize',15);

legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','Discrete Obs on Bdy','Discrete Obs on Bdy Diag Jumps','location','southeast');
xlabel('log_2(1 / mean path length)')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;
ca.XLim(1) = 1;

%% Publication: Continuous vs discrete vs discrete w/ diag corrected
saveOn = 0;

% continuous data
load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
results_cont = results;
sim_spec = [results_cont.sim_spec];

pathLen_cont = [sim_spec.path_len];
CI = [results_cont.Deff95CI];
CIlow = CI(1:2:end);
CIhigh = CI(2:2:end);
Deff_cont = [results_cont.Deff];

% discrete data
load('Results/results_2d_square_2017_07_19_09_37_30.mat','results_homog');
results_discrete = results_homog;

m = 2.^[2:8];
pathLen_disc = 1./m;%1./[ghParams_discrete.m];
Deff_disc = [results_discrete.Deff];

% discrete diag data. results collected at two points in time. methodologies for diag jumps
% differed.
load('Results_2d_square_diagJumpsCorrected/results_2017_08_13_18_38_51.mat','results_homog');
results_discrete_diagC = results_homog(1:8);
Deff_disc_diagC = [results_discrete_diagC.Deff];

% newer discrete diag jump data. methodology differed slightly
%{
load('Results/results_2d_square_diagJumps_2017_05_20_12_22_53.mat','results_homog');
results_discrete_diagC = results_homog;
results_discrete_diagC = results_discrete_diagC(1:16);
Deff_disc_diagC = [results_discrete_diagC.Deff];
Deff_disc_diagC = Deff_disc_diagC(2:2:end)./Deff_disc_diagC(1:2:end);
%}
m_C = 2.^[2:9];
pathLen_disc_diagC = (.5*(sqrt(2) + 1))./m_C;%(.5*(sqrt(2) + 1))./[ghParams_discrete_diagC.m];

% plot
fh = figure;
hold on
errorbar(log2(1./pathLen_cont),Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',30);
h = plot(log2(1./pathLen_disc),Deff_disc,'s-','markersize',10);
set(h, 'MarkerFaceColor', get(h, 'Color'));
h = plot(log2(1./pathLen_disc_diagC),Deff_disc_diagC,'d-','markersize',10);
set(h, 'MarkerFaceColor', get(h, 'Color'));

legend('Continuous','Discrete','Discrete w/ diagonal jumps','location','southeast');
xlabel('Mean path length')
ylabel('D_e')
ca = gca;
ca.FontSize = 18;
ca.XLim(1) = 1;
ca.XTickLabel = {'2^{-2}','2^{-3}','2^{-4}','2^{-5}','2^{-6}','2^{-7}','2^{-8}'};
axis([1.5 8 .66 .78])
if saveOn

    dirName = '';
    filename = 'contVsDiscVsDiscdiag';
    
    mySaveFig([dirName filename],fh,'fig')
    mySaveFig([dirName filename],fh,'png')
    mySaveFig([dirName filename],fh,'eps')
    
end