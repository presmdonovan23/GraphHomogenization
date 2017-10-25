% runs square obstruction with

saveOn = 0;
plotOn = 1;
rho = .5;
mVals = 2.^[2];%2.^[2:9];
geometryName = 'squareBdySlow';
delta = 2;
dim = 2;
D0 = 1;

numTraj = 0;
startNodeInd = [];

rate = [];

%% no diag jumps
diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps

i = 1;
for m = mVals
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    rateCoeffs.alpha = 0;
    rateCoeffs.K1 = [];
    rateCoeffs.K2 = [];
    rateCoeffs.delta = delta;
    
    ctr = .75 - .5*h;
    %ctr = .5 - .5*h;
    rhoCorrected = rho - .5*h;

    geometry.dim = dim;
    geometry.name = geometryName;
    geometry.m = m;
    geometry.rho = rhoCorrected;
    geometry.diagJumps = diagJumps;
    geometry.ctr = ctr;
    geometry.specialSetting_m2 = [];
    geometry.rateCoeffs = rateCoeffs;

    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rhoCorrected,m,rate,rateCoeffs,diagJumps,ctr);
    results(i) = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry);
    
    if plotOn
        plotEdges = 1;
        plotObs = 1;
        plotRates = 0;
        
        obRad = rho/2;
        
        plotGraph(nodes,edges,edgeRates,edgeJumps,geometry.name,obRad,geometry.ctr,plotEdges,plotObs,plotRates);
    end
    i = i+1;
end

if saveOn
    saveResults(results);
end

%% diag jumps
diagJumps = 1; % diagJumps = 2 => corrected diagonal jumps

i = 1;
for m = mVals
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    rateCoeffs.alpha = 0;
    rateCoeffs.K1 = [];
    rateCoeffs.K2 = [];
    rateCoeffs.delta = delta;
    
    ctr = .75 - .5*h;
    %ctr = .5 - .5*h;
    rhoCorrected = rho - .5*h;

    geometry.dim = dim;
    geometry.name = geometryName;
    geometry.m = m;
    geometry.rho = rhoCorrected;
    geometry.diagJumps = diagJumps;
    geometry.ctr = ctr;
    geometry.specialSetting_m2 = [];
    geometry.rateCoeffs = rateCoeffs;

    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rhoCorrected,m,rate,rateCoeffs,diagJumps,ctr);
    results(i) = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry);
    
    if plotOn
        plotEdges = 1;
        plotObs = 1;
        plotRates = 0;
        
        obRad = rho/2;
        
        plotGraph(nodes,edges,edgeRates,edgeJumps,geometry.name,obRad,geometry.ctr,plotEdges,plotObs,plotRates);
    end
    i = i+1;
end

if saveOn
    saveResults(results);
end
