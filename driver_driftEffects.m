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
plotOn = 1;
rate = []; % field is current obsolete
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

%% regular
geometry.name = 'square';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,obRad,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{1} = saveResults(results);
end

%% bonding
rateCoeffs.delta = .5;
geometry.name = 'squareBonding';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn

    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,obRad,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{2} = saveResults(results);
end

%% boundary attraction
rateCoeffs.delta = 2;
geometry.name = 'squareBdyAttract';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,obRad,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{3} = saveResults(results);
end

%% boundary repel
rateCoeffs.delta = 2;
geometry.name = 'squareBdyRepel';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,obRad,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{4} = saveResults(results);
end