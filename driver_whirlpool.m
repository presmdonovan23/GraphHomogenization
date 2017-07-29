removeEdges = 0;
rateReduction = .5;
numTraj = 00;
plotOn = 1;
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

[ghInput, results_homog, results_mc] = getDeff_whirlpool(ghParams,numTraj,removeEdges,rateReduction,plotOn);