% Plot results from square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary
saveOn = 1;
%% Continuous vs discrete vs discrete w/ diag vs discrete w/ diag corrected
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
ca.XTickLabel = {'Simple','Bonding','Repulsion','Attraction'};
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