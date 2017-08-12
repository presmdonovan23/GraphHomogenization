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
rate = []; % field is current obsolete
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);

[ghInput, results_homog, results_mc] = getDeff_whirlpool(ghParams,numTraj,removeEdges,rateReduction,plotOn,plotLambda,plotNodeLabels);

dirname = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';

fh = gcf;
mySaveFig([dirname 'detailedBalance2'],fh,'fig')
mySaveFig([dirname 'detailedBalance2'],fh,'png')
mySaveFig([dirname 'detailedBalance2'],fh,'eps')