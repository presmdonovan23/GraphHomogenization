% runs square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary

plotOn = 1;
saveOn = 0;

% geometry
dim = 2;
m = 8;
h = 1/m;
name = 'square';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 0;
driftMult = 0;
driftDecay = [];

% neutral
specialSetting = 'none';
obSlowdownFctr = [];
bdyDist = [];

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
            
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

if saveOn
    saveresults(results);
end
% boundary bonding
specialSetting = 'bdyBonding';
obSlowdownFctr = 1/2;
bdyDist = h;

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
            
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

if saveOn
    saveresults(results);
end
% boundary attraction
specialSetting = 'bdyAttractRepel';
obSlowdownFctr = 1/2;
bdyDist = [];

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
            
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

if saveOn
    saveresults(results);
end
% boundary repulsion
specialSetting = 'bdyAttractRepel';
obSlowdownFctr = 2;
bdyDist = [];

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
            
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

if saveOn
    saveresults(results);
end