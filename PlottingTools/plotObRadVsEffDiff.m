function fh = plotObRadVsEffDiff(results,fh)

if nargin < 2 || isempty(fh)
    fh = figure;
end

figure(fh);
nResults = length(results);

geometry = [results.geometry];
mc = [results.mc];

hold on

plot([geometry.obRad],[results.Deff],'.-','markersize',15)

if mc(1).numTraj > 0
    CI = reshape([mc.Deff_95CI],2,nResults)';
    errorbar([geometry.obRad],[mc.Deff],...
        .5*(CI(:,2) - CI(:,1)),'.-','markersize',15);
end

xlabel('Obstruction Radius')
ylabel('D_e')

legend('Homogenization Theory','Monte Carlo (95% CI)')

end