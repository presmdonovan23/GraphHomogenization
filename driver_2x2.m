numTraj = 5000;
setting = 'blockOneSite';
%setting = 'slowOneSite';
%setting = '';
delta = .1; % only required if setting = slowOneSite

% DO NOT CHANGE
rho = 0;
m = 2;
dim = 2;
geometry = 'square';
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

[ghInput, results_homog, results_mc] = getDeff_2x2(ghParams,numTraj,setting,delta);