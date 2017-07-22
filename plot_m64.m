%% uncomment lines below to include monte carlo

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

mysavefig('Deff_m64_drift',fh,'png')
mysavefig('Deff_m64_drift',fh,'fig')
mysavefig('Deff_m64_drift',fh,'eps')