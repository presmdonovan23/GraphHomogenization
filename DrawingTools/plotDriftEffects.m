%% interaction only at boundary
% Plot results from square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary
saveOn = 1;

dirName = '/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/';

dirnameSub = 'Results_2d_square/';
filename = 'results_2017_08_13_10_24_53.mat';
load([dirName dirnameSub filename]);
results_regular = results_homog;
ghParams_regular = ghParams;

dirnameSub = 'Results_2d_squareBonding/';
filename = 'results_2017_08_13_10_24_53.mat';
load([dirName dirnameSub filename]);
results_bonding = results_homog;
ghParams_bonding = ghParams;

dirnameSub = 'Results_2d_squareBdyRepel/';
filename = 'results_2017_08_13_10_24_54.mat';
load([dirName dirnameSub filename]);
results_bdyRepel = results_homog;
ghParams_bdyRepel = ghParams;

dirnameSub = 'Results_2d_squareBdyAttract/';
filename = 'results_2017_08_13_10_24_54.mat';
load([dirName dirnameSub filename]);
results_bdyAttract = results_homog;
ghParams_bdyAttract = ghParams;

fh = figure;
bh = bar([results_regular.Deff,results_bonding.Deff,results_bdyRepel.Deff,results_bdyAttract.Deff]);
ca = gca;
ca.XTickLabel = {'Neutral','Bonding','Repulsion','Attraction'};
ylabel('D_e')
ca = gca;
ca.FontSize = 18;

text(bh.XData,bh.YData',num2str(bh.YData','%0.4f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'fontsize',18)

if saveOn
    dirName = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';
    filename = 'driftEffectsDeff';
    mySaveFig([dirName 'driftEffectsDeff'],fh,'fig')
    mySaveFig([dirName 'driftEffectsDeff'],fh,'png')
    mySaveFig([dirName 'driftEffectsDeff'],fh,'eps')
end

%% interaction via drift field
% neutral vs attraction vs repulsion (drift field)
saveOn = 0;

load('Results/results_2d_circle_m64_alpha0_noapprox.mat')
% need to run monte carlo for this case
%load('results_2d_circle_m64_alpha0_noapprox_10000.mat')
results_alpha0 = results;
load('Results/results_2d_circle_m64_alpha1_noapprox_10000.mat')
results_alpha1 = results;
load('Results/results_2d_circle_m64_alpham1_noapprox_10000.mat')
results_alpham1 = results;

Deff_alpha0 = [1; zeros(10,1)];
Deff_alpha0_mc = [1; zeros(10,1)];
Deff_alpha0_mc_ci = zeros(11,2);

Deff_alpha1 = [1; zeros(10,1)];
Deff_alpha1_mc = [1; zeros(10,1)];
Deff_alpha1_mc_ci = zeros(11,2);

Deff_alpham1 = [1; zeros(10,1)];
Deff_alpham1_mc = [1; zeros(10,1)];
Deff_alpham1_mc_ci = zeros(11,2);

for i = 2:10
    idx = i-1;
    Deff_alpha0(i) = results_alpha0(idx).Deffres.Deff_homog(1);
    Deff_alpha0_mc(i) = results_alpha0(idx).Deffres.Deff_MC(1);
    Deff_alpha0_mc_ci(i,:) = results_alpha0(idx).Deffres.Deff_MC_95CI;

    Deff_alpha1(i) = results_alpha1(idx).Deffres.Deff_homog(1);
    Deff_alpha1_mc(i) = results_alpha1(idx).Deffres.Deff_MC(1);
    Deff_alpha1_mc_ci(i,:) = results_alpha1(idx).Deffres.Deff_MC_95CI;

    Deff_alpham1(i) = results_alpham1(idx).Deffres.Deff_homog(1);
    Deff_alpham1_mc(i) = results_alpham1(idx).Deffres.Deff_MC(1);
    Deff_alpham1_mc_ci(i,:) = results_alpham1(idx).Deffres.Deff_MC_95CI;
end

markersize = 20;
rho = (0:.1:1);
R = rho./2;
fh = figure;
hold on
plot(R,Deff_alpha0,'k.-','markersize',markersize)
plot(R,Deff_alpha1,'b.-','markersize',markersize)
%errorbar(R,Deff_alpha1_mc,Deff_alpha1_mc_ci(:,2) - Deff_alpha1_mc_ci(:,1),...
%        'sb','markersize',markersize/3);
plot(R,Deff_alpham1,'r.-','markersize',markersize)
%errorbar(R,Deff_alpham1_mc,Deff_alpham1_mc_ci(:,2) - Deff_alpham1_mc_ci(:,1),...
%        'sr','markersize',markersize/3);
xlabel('R')
ylabel('D_{eff}')
%legend('No force','Attraction (homogenization)','Attraction (Monte Carlo)',...
%    'Repulsion (homogenization)','Repulsion (Monte Carlo)',...
%    'location','southwest')
legend( 'No force','Attraction (homogenization)','Repulsion (homogenization)',...
        'location','southwest')


ca = gca;
ca.FontSize = 18;

if saveOn
    
    mysavefig('Deff_m64_drift',fh,'png')
    mysavefig('Deff_m64_drift',fh,'fig')
    mysavefig('Deff_m64_drift',fh,'eps')
    
end