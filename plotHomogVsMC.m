function fh = plotHomogVsMC(results_homog,results_mc)

fh = figure;
hold on
CI = reshape([results_mc.Deff_95CI],2,3)';
plot([ghParams.rho],[results_homog.Deff],'.-','markersize',15)
errorbar([ghParams.rho],[results_mc.Deff],...
    .5*(CI(:,2) - CI(:,1)),'.-','markersize',15);

xlabel('rho')
ylabel('Deff')

legend('Homogenization Theory','Monte Carlo (95% CI)')

end