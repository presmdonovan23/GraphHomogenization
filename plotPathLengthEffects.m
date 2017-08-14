%% Continuous homogenization vs discrete homogenization
load('Results_2d_square/results_2017_08_05_10_35_27.mat')

x_homog = [1:9];
y_homog = [2/3 results_homog.Deff];
x_cont = linspace(1,9,100);
y_cont = .77*ones(1,100);

fh = figure;

hold on
plot(x_cont,y_cont,'r-','linewidth',2)
plot(x_homog,y_homog,'b.-','markersize',20,'linewidth',2)

axis([1 9 .66 .78]);
axis square

xlabel('Path Length (h)');
ylabel('D_e^h');

lh = legend('Continuous Homogenization','Graph Homogenization','location','southeast');

ca = gca;
ca.FontSize = 18;
ca.XTick = [1:9];
for i = 1:9
    ca.XTickLabel{i} = sprintf('2^{-%d}',i);
end

dirname = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';

%mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'fig')
%mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'png')
%mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'eps')

%% Plot graph setting w/ square in upper right corner (m = 2,4,8,16)
fh = figure;
i = 1;
ctr = .75;
saveOn = 0;
drawObs = 1;
for m = [2 4 8 16]
    sph(i) = subplot(2,2,i);
    if i == 2 || i == 4
        sph(i).Position = sph(i).Position - [.15 .1 -.15 -.15];
    else
        sph(i).Position = sph(i).Position - [.1 .1 -.15 -.15];
    end
    rho = .5;
    %m = 2;
    geometry = 'square';
    drawObs = 1;
    drawSetting(rho,m,ctr,geometry,drawObs,saveOn,fh);
    i=i+1;
end

dirname = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';

%mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'fig')
%mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'png')
%mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'eps')

%% Continuous vs discrete vs discrete w/ diag vs discrete w/ diag corrected
dirname = '/Users/prestondonovan/Documents/School/Research/MATLAB Code/Discrete Homogenization/GraphHomogenization/';
load('/Users/prestondonovan/Documents/School/Research/MATLAB Code/Data Processing/Simulation_2d/FixedPathLength/results.mat');
results_cont = results;

load([dirname 'results/results_2d_square_2017_07_19_09_37_30.mat']);
results_discrete = results_homog;
ghParams_discrete = ghParams;
%ghInput_discrete = ghInput;

load([dirname 'Results_2d_square_diagJumpsCorrected/results_2017_08_13_18_38_51.mat']);
results_discrete_diagC = results_homog;
ghParams_discrete_diagC = ghParams;
ghInput_discrete_diagC = ghInput;

load([dirname 'Results_2d_square_diagJumps/results_2017_07_29_13_28_21.mat']);
results_discrete_diag = results_homog;
ghParams_discrete_diag = ghParams;
ghInput_discrete_diag = ghInput;

load([dirname 'Results_2d_squareBdySlow/results_2017_08_14_12_49_25.mat'])
results_obsBdy = results_homog;
ghInput_obsBdy = ghInput;
ghParams_obsBdy = ghParams;

load([dirname 'Results_2d_squareBdySlow_diagJumps/results_2017_08_14_12_49_38.mat'])
results_obsBdy_diag = results_homog;
ghInput_obsBdy_diag = ghInput;
ghParams_obsBdy_diag = ghParams;

sim_spec = [results_cont.sim_spec];

pathLen_cont = [sim_spec.path_len];
pathLen_disc = 1./[ghParams_discrete.m];
pathLen_disc_diagC = (.5*(sqrt(2) + 1))./[ghParams_discrete_diagC.m];
pathLen_disc_diag = (.5*(sqrt(2) + 1))./[ghParams_discrete_diag.m];
pathLen_obsBdy = 1./[ghParams_obsBdy.m];
pathLen_obsBdy_diag = (.5*(sqrt(2) + 1))./[ghParams_obsBdy_diag.m];

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
%{
subplot(1,3,1)
hold on
errorbar(pathLen_cont,Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(pathLen_disc,Deff_disc,'.-','markersize',15);
plot(pathLen_disc_diag,Deff_disc_diag,'.-','markersize',15);
plot(pathLen_disc_diagC,Deff_disc_diagC,'.-','markersize',15);
plot(pathLen_obsBdy,Deff_obsBdy,'.-','markersize',15);
plot(pathLen_obsBdy_diag,Deff_obsBdy_diag,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','Discrete Obs on Bdy','Discrete Obs on Bdy Diag Jumps','location','southeast');
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
plot(1./pathLen_obsBdy,Deff_obsBdy,'.-','markersize',15);
plot(1./pathLen_obsBdy_diag,Deff_obsBdy_diag,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','Discrete Obs on Bdy','Discrete Obs on Bdy Diag Jumps','location','southeast');
xlabel('1 / mean path length')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;
%}
%subplot(1,3,3)
hold on
errorbar(log2(1./pathLen_cont),Deff_cont,...
    .5*(CIhigh - CIlow),'.-','markersize',15);
plot(log2(1./pathLen_disc),Deff_disc,'.-','markersize',15);
plot(log2(1./pathLen_disc_diag),Deff_disc_diag,'.-','markersize',15);
plot(log2(1./pathLen_disc_diagC),Deff_disc_diagC,'.-','markersize',15);
plot(log2(1./pathLen_obsBdy),Deff_obsBdy,'.-','markersize',15);
plot(log2(1./pathLen_obsBdy_diag),Deff_obsBdy_diag,'.-','markersize',15);

lh = legend('Continuous','Discrete','Discrete Diag Jumps','Discrete Diag Jumps (corrected)','Discrete Obs on Bdy','Discrete Obs on Bdy Diag Jumps','location','southeast');
xlabel('log(1 / mean path length)')
ylabel('Deff')
ca = gca;
ca.FontSize = 18;
ca.XLim(1) = 1;