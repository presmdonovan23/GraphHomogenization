%% A simple driver
plotOn = 1;
saveOn = 0;

% geometry parameters
dim = 2;
m = 5;
name = 'circle';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

% monte carlo parameters
startNodeInd = 1;
numTraj = 100;

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);
                
if saveOn
    saveResults(results);
end