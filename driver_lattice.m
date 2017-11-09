% misc
saveOn = 0;
plotOn = 1;

% geometry
mVals = 2.^[3];
obRadVals = [.125 .25];

dim = 2;
name = 'square';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

% monte carlo
startNodeInd = 1;
numTraj = 100;

for i = 1:length(mVals)
    m = mVals(i);
    for j = 1:length(obRadVals)
        obRad = obRadVals(j);
        
        latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                    obCtr, diagJumps, specialSetting, ...
                    driftMult, driftDecay, obSlowdownFctr, bdyDist);

        [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

        curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

        results(i,j) = curRes;
    end
end

if plotOn
    for i = 1:length(mVals)
    
        plotObRadVsEffDiff(results(i,:));
        
    end
end

if saveOn
    saveResults(results);
end