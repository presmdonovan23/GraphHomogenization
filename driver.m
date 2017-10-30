%% A simple driver
plotOn = 1;

% geometry parameters
dim = 2;
geometryName = 'circle';
obRad = .25;
m = 5;
diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
obCtr = .5;
specialSetting_m2 = 'none';

% monte carlo parameters
startNodeInd = 1;
numTraj = 100;

rate = []; % field is current obsolete

% rate parameters
alpha = 0;
K1 = 25;
K2 = 10;
delta = 1.5;
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;
rateCoeffs.delta = delta;

geometry.dim = dim;
geometry.name = geometryName;
geometry.m = m;
geometry.obRad = obRad;
geometry.diagJumps = diagJumps;
geometry.ctr = obCtr;
geometry.specialSetting_m2 = specialSetting_m2;
geometry.rateCoeffs = rateCoeffs;
        
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(...
    dim,geometryName,obRad,m);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,...
                    numTraj,startNodeInd,geometry);