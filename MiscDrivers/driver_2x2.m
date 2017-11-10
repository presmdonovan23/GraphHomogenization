%% blockOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = .25;
obCtr = [.75 .75];
diagJumps = 0;
specialSetting = [];
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

% Monte Carlo
numTraj = 10000;
startNodeInd = [];

% get Deff
latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

%% slowOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = .25;
obCtr = [.75 .75];
diagJumps = 0;
specialSetting = 'm2_slowOneSite'; %slowOneSite
driftMult = 0;
driftDecay = [];
obSlowdownFctr = .5;
bdyDist = [];

% Monte Carlo
numTraj = 10000;
startNodeInd = [];

% get Deff
latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);