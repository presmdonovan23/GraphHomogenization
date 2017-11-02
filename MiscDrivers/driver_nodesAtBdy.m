%% bdySlow (square)

saveOn = 0;

% geometry
mVals = 2.^[2:9];

dim = 2;
name = 'square';
obRad = .25;
diagJumps = 0;
specialSetting = 'bdySlow';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = 2;

% monte carlo
startNodeInd = 1;
numTraj = 100;

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    obCtr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    bdyDist = .9999*h;
    
    latticeGeo = LatticeGeometry(dim, m, name, obRadCorrected, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

    results(i) = curRes;

end

if saveOn
    saveResults(results);
end

%% bdySlow with diagonal jumps (square)

saveOn = 0;

% geometry
mVals = 2.^[2:9];

dim = 2;
name = 'square';
obRad = .25;
diagJumps = 1;
specialSetting = 'bdySlow';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = 2;

% monte carlo
startNodeInd = 1;
numTraj = 100;

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    obCtr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    bdyDist = .9999*h;
    
    latticeGeo = LatticeGeometry(dim, m, name, obRadCorrected, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

    results(i) = curRes;

end

if saveOn
    saveResults(results);
end