
saveOn = 0;

% geometry parameters
dim = 2;
geometryName = 'square';
D0 = 1;
rhoVals = [.25 .5 .75];
mVals = 2.^[3];
diagJumps = 2; % diagJumps = 2 => corrected diagonal jumps
ctr = .5;
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

idx = 1;
for m = mVals
    for rho = rhoVals
        
        geometry.dim = dim;
        geometry.name = geometryName;
        geometry.m = m;
        geometry.rho = rho;
        geometry.diagJumps = diagJumps;
        geometry.ctr = ctr;
        geometry.specialSetting_m2 = specialSetting_m2;
        geometry.rateCoeffs = rateCoeffs;
        
        [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);

        curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry);
        
        results(idx) = curRes;
        idx = idx+1;
    
    end
end

if saveOn
    saveResults(results);
end