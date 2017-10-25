
saveOn = 0;
plotOn = 1;

% geometry parameters
dim = 2;
geometryName = 'square';
D0 = 1;
rhoVals = [.25 .5 .75];
mVals = 2.^[3 4];
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

clear results;
for i = 1:length(mVals)
    m = mVals(i);
    for j = 1:length(rhoVals)
        rho = rhoVals(j);
        
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
        
        results(i,j) = curRes;
        
    end
    
end

if plotOn
    for i = 1:length(mVals)
    
        plotRhoVsEffDiff(results(i,:));
        
    end
end

if saveOn
    saveResults(results);
end