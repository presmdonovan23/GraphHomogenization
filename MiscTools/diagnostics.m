%%
% The purpose of this script is to ensure that, after doing some
% development, the output has not changed. We check the results against
% some files that serve as the gold standard for correctness. These files
% will be included in a future release.
%% (circle)
mVals = 2.^[2:8];

dim = 2;
name = 'circle';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(mVals)
    m = mVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

load('Results/results_2d_circle_2017_07_19_09_30_19.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

%% (3D circle)

obRadVals = [0  .35 .65];

dim = 3;
name = 'circle';
m = 16;
obCtr = [.5 .5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(obRadVals)
    obRad = obRadVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

load('Results/results_3d_circle_2017_06_03_15_13_54.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

%% (3D square)

obRadVals = [.05:.05:.45];

dim = 3;
name = 'square';
m = 16;
obCtr = [.5 .5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(obRadVals)
    obRad = obRadVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

load('Results/results_3d_square_2017_06_03_12_14_50.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

%% bdyBonding, bdyAttract, bdyRepel, (square)
% runs square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary

dim = 2;
m = 8;
name = 'square';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);


% regular

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

abs(results.Deff - .743103104943094)

% bonding
latticeGeo.specialSetting = 'bdyBonding';
latticeGeo.obSlowdownFctr = 1/2;
latticeGeo.bdyDist = 1/m;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

abs(results.Deff - .524543368195125)

% boundary attraction
latticeGeo.specialSetting = 'bdyAttractRepel';
latticeGeo.obSlowdownFctr = 1/2;
latticeGeo.bdyDist = [];

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

abs(results.Deff - .705490976286049)

% boundary repel
latticeGeo.specialSetting = 'bdyAttractRepel';
latticeGeo.obSlowdownFctr = 2;
latticeGeo.bdyDist = [];

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

abs(results.Deff - .680616147421872)

%% none (square)

mVals = 2.^[2:8];

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

for i = 1:length(mVals)
    m = mVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

load('Results/results_2d_square_2017_07_19_09_37_30.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])
%% diagonal jumps corrected (square)

mVals = 2.^[2:9];
dim = 2;
name = 'square';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 2;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(mVals)
    m = mVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

% to approximate this data, see comment in rate_lattice.m (i.e., make jump
% rates along border scaled by 1.5 rather than 2)
load('Results_2d_square_diagJumpsCorrected/results_2017_08_13_18_38_51.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

load('Results_2d_square_diagJumpsCorrected/results_2017_07_31_19_35_21.mat','results_homog');
norm([results(1:7).Deff] - [results_homog.Deff])

%% diagonal jumps (square)

mVals = 2.^[2:8];

dim = 2;
name = 'square';
obRad = .25;
obCtr = [.5 .5];
diagJumps = 1;
specialSetting = 'none';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(mVals)
    m = mVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;
       
end

load('Results_2d_square_diagJumps/results_2017_07_29_13_28_21.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

%% bdySlow (square)

mVals = 2.^[2:9];

dim = 2;
name = 'square';
obRad = .25;
diagJumps = 0;
specialSetting = 'bdySlow';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    obCtr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    
    
    latticeGeo = LatticeGeometry(dim, m, name, obRadCorrected, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;

end

load('Results_2d_squareBdySlow/results_2017_08_14_12_49_25.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])

%% bdySlow with diagonal jumps (square)

mVals = 2.^[2:9];
dim = 2;
name = 'square';
obRad = .25;
diagJumps = 1;
specialSetting = 'bdySlow';
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(mVals)
    m = mVals(i);
    h = 1/m;
    
    obCtr = .75 - .5*h;
    obRadCorrected = obRad - .25*h;
    
    latticeGeo = LatticeGeometry(dim, m, name, obRadCorrected, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;

end

load('Results_2d_squareBdySlow_diagJumps/results_2017_08_14_12_49_38.mat','results_homog');
norm([results.Deff] - [results_homog.Deff])
%% drift field repel (circle)

m = 64;
dim = 2;
name = 'circle';
obRadVals = (.05:.05:.45);
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = 25;
driftDecay = 10;
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(obRadVals)
    obRad = obRadVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;

end

results_orig = load('Results/results_2d_circle_m64_alpha1_noapprox.mat');
results_orig = [results_orig.results.Deffres];
Deff_orig = [results_orig.Deff_homog];
Deff_orig = Deff_orig(1,1:2:end);
norm([results.Deff] - Deff_orig)

%% drift field attract (circle)

m = 64;
dim = 2;
name = 'circle';
obRadVals = (.05:.05:.45);
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'none';
driftMult = -25;
driftDecay = 10;
obSlowdownFctr = [];
bdyDist = [];

for i = 1:length(obRadVals)
    obRad = obRadVals(i);
    
    latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
                obCtr, diagJumps, specialSetting, ...
                driftMult, driftDecay, obSlowdownFctr, bdyDist);
    
    [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

    curRes = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],latticeGeo);

    results(i) = curRes;

end

results_orig = load('Results/results_2d_circle_m64_alpham1_noapprox.mat');
results_orig = [results_orig.results.Deffres];
Deff_orig = [results_orig.Deff_homog];
Deff_orig = Deff_orig(1,1:2:end);
norm([results.Deff] - Deff_orig)

%% blockOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = .25;
obCtr = [.75 .75];
diagJumps = 0;
specialSetting = [];
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

% Monte Carlo
numTraj = 10000;
startNodeInd = [];

% get Deff
latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

norm(results.Deff - 2/3)
norm(results.mc.Deff - 2/3)
%% slowOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = .25;
obCtr = [.75 .75];
diagJumps = 0;
specialSetting = 'm2_slowOneSite'; %slowOneSite
driftMult = 0;
driftDecay = [];
obSlowdownFctr = .35;
bdyDist = [];

% Monte Carlo
numTraj = 10000;
startNodeInd = [];

% get Deff
latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo);

correctDeff = (.5 + .5*obSlowdownFctr)
norm(results.Deff - correctDeff)
norm(results.mc.Deff - correctDeff)