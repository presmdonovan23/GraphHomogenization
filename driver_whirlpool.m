removeEdges = 1;
rateReduction = .5;
numTraj = 00;
plotOn = 1;
plotLambda = 0;
plotNodeLabels = 1;
% DO NOT CHANGE
rho = .25;
m = 3;
dim = 2;
geometry = 'circle';
diagJumps = 0;
D0 = 1;
alpha = 0;
K1 = 25;
K2 = 10;
rate = [];
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);

ghInput = GraphHomogInput(ghParams);
ghInput = getDeff_whirlpool(ghInput,removeEdges,rateReduction,plotOn,plotLambda,plotNodeLabels);

results_homog = getDeff_homog(ghInput);

for i = 1:8
    startNodeInd(1 + (i-1)*numTraj:i*numTraj) = i;
end
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );

dirname = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';

if plotOn
    
    plotEdges = 1;
    plotObs = 0;
    fh = drawWhirlpoolSetting( ghParams, ghInput, plotEdges, plotObs );

    %mySaveFig([dirname 'detailedBalance2'],fh,'fig')
    %mySaveFig([dirname 'detailedBalance2'],fh,'png')
    %mySaveFig([dirname 'detailedBalance2'],fh,'eps')
    
end
