%% neutral vs bonding vs attraction vs repulsion
% runs square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary

saveOn = 0;
obRad = .25;
m = 8;
h = 1/m;
dim = 2;

diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
ctr = .5;
specialSetting_m2 = 'none';

startNodeInd = 1;
numTraj = 0;
plotOn = 0;
drawObs = 1;
drawEdges = 1;
drawRates = 1;

rateCoeffs.dist = .9999*h;
rateCoeffs.alpha = 0;
rateCoeffs.K1 = [];
rateCoeffs.K2 = [];
geometry.dim = dim;
geometry.m = m;
geometry.obRad = obRad;
geometry.diagJumps = diagJumps;
geometry.ctr = ctr;
geometry.specialSetting_m2 = specialSetting_m2;

% regular
geometry.name = 'square';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

abs(results.Deff - .7431)

% bonding
rateCoeffs.delta = .5;
geometry.name = 'squareBonding';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

abs(results.Deff - .5245)

% boundary attraction
rateCoeffs.delta = 2;
geometry.name = 'squareBdyAttract';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

abs(results.Deff - .7055)

% boundary repel
rateCoeffs.delta = 2;
geometry.name = 'squareBdyRepel';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

abs(results.Deff - .6806)

%% discrete


dim = 2;
obRad = .25;
ctr = [];
specialSetting_m2 = [];

load('Results/results_2d_square_2017_07_19_09_37_30.mat','results_homog');
geometryName = 'square';
mVals = 2.^[2:8];
diagJumps = 0;

rateCoeffs.alpha = 0;
rateCoeffs.K1 = 0;
rateCoeffs.K2 = 0;
rateCoeffs.delta = 0;

for i = 1:length(mVals)
    m = mVals(i);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],[]);

    results(i) = curRes;
       
end

norm([results.Deff] - [results_homog.Deff])
%% discrete with diag corrected

load('Results_2d_square_diagJumpsCorrected/results_2017_08_13_18_38_51.mat','results_homog');

%% discrete with diag


dim = 2;
obRad = .25;
ctr = [];
specialSetting_m2 = [];

load('Results_2d_square_diagJumps/results_2017_07_29_13_28_21.mat','results_homog');
geometryName = 'square';
mVals = 2.^[2:8];
diagJumps = 1;

rateCoeffs.alpha = 0;
rateCoeffs.K1 = 0;
rateCoeffs.K2 = 0;
rateCoeffs.delta = 0;

for i = 1:length(mVals)
    m = mVals(i);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],[]);

    results(i) = curRes;
       
end

norm([results.Deff] - [results_homog.Deff])

%% slowed rates at boundary

dim = 2;
obRad = .25;
specialSetting_m2 = [];

load('Results_2d_squareBdySlow/results_2017_08_14_12_49_25.mat','results_homog');
geometryName = 'squareBdySlow';
mVals = 2.^[2:9];
diagJumps = 0;

rateCoeffs.alpha = 0;
rateCoeffs.K1 = 0;
rateCoeffs.K2 = 0;
rateCoeffs.delta = 2;

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    
    ctr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,obRadCorrected,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],[]);

    results(i) = curRes;

end

norm([results.Deff] - [results_homog.Deff])

%% slowed rates at boundary with diagonal jumps


dim = 2;
obRad = .25;
specialSetting_m2 = [];

load('Results_2d_squareBdySlow_diagJumps/results_2017_08_14_12_49_38.mat','results_homog');
geometryName = 'squareBdySlow';
mVals = 2.^[2:9];
diagJumps = 1;

rateCoeffs.alpha = 0;
rateCoeffs.K1 = 0;
rateCoeffs.K2 = 0;
rateCoeffs.delta = 2;

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    
    ctr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometryName,obRadCorrected,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],[]);

    results(i) = curRes;

end

norm([results.Deff] - [results_homog.Deff])
%% drift field

load('Results/results_2d_circle_m64_alpha1_noapprox.mat')
results_orig = [results.Deffres];
Deff_orig = [results_orig.Deff_homog];
Deff_orig = Deff_orig(1,1:2:end);

clear results

saveOn = 0;
obRads = (.05:.05:.45);
m = 64;
h = 1/m;
dim = 2;

diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
ctr = .5;
specialSetting_m2 = 'none';

startNodeInd = 1;
numTraj = 0;
plotOn = 0;
drawObs = 1;
drawEdges = 1;
drawRates = 1;

rateCoeffs.dist = .9999*h;
rateCoeffs.alpha = 1;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;
geometry.dim = dim;
geometry.m = m;

geometry.diagJumps = diagJumps;
geometry.ctr = ctr;
geometry.specialSetting_m2 = specialSetting_m2;

geometry.name = 'circle';
geometry.rateCoeffs = rateCoeffs;

for i = 1:length(obRads)
    obRad = obRads(i);
    geometry.obRad = obRad;
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rateCoeffs,diagJumps,ctr,specialSetting_m2);
    
    results(i) = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);
    
end

norm([results.Deff] - Deff_orig)