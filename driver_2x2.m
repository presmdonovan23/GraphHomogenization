numTraj = 0;
specialSetting_m2 = 'blockOneSite';
specialSetting_m2 = 'slowOneSite';
delta = .1; % only required if setting = slowOneSite

% DO NOT CHANGE
rho = 0;
m = 2;
dim = 2;
geometry = 'square';
diagJumps = 0;
D0 = 1;
rate = [];
rateCoeffs.alpha = 0;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;
rateCoeffs.delta = delta;
ctr = [];
negateSteps = 1;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
ghInput = GraphHomogInput(ghParams);
results_homog = getDeff_homog_D2_m2(ghInput, specialSetting_m2);
results_MC = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn, negateSteps );