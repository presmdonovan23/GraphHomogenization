%% square drift test

% geometry parameters
dim = 2;
geometryName = 'square';
obRadVals = [.125 .25 .375];
mVals = 2.^[3 4];
diagJumps = 2; % diagJumps = 2 => corrected diagonal jumps
ctr = .5;
specialSetting_m2 = 'none';

% monte carlo parameters
startNodeInd = 1;
numTraj = 100;

rate = []; % field is current obsolete

% rate parameters
rateCoeffs.alpha = 1;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;
rateCoeffs.delta = 1.5;

clear results;
for i = 1:length(mVals)
    m = mVals(i);
    for j = 1:length(obRadVals)
        obRad = obRadVals(j);
        
        geometry.dim = dim;
        geometry.name = geometryName;
        geometry.m = m;
        geometry.obRad = obRad;
        geometry.diagJumps = diagJumps;
        geometry.ctr = ctr;
        geometry.specialSetting_m2 = specialSetting_m2;
        geometry.rateCoeffs = rateCoeffs;
        
        [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);

        curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry);
        
        results(i,j) = curRes;
        
    end
    
end
