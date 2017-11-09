%% blockOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = 0;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'm2_blockOneSite'; %slowOneSite
driftMult = 0;
driftDecay = [];
obSlowdownFctr = [];
bdyDist = [];

% Monte Carlo
numTraj = 0000;
startNodeInd = [];
negateSteps = 1;
plotOn = 0;

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

pi0 = (1/3)*ones(3,1);
pi0_relRes = [];
pi0_flag = [];
pi0_time = [];

unitCell_soln = pi0;
unitCell_RHS = [];
unitCell_relRes = [];
unitCell_solvability = [];
unitCell_flag = [];
unitCell_time = [];

LUtime = [];

[Deff, DeffTerm1, DeffTerm2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln );

results.L =                     L;
results.nodes =                 nodes;
results.edges =                 edges;
results.edgeRates =             edgeRates;
results.edgeJumps =             edgeJumps;
results.LU_time =               LUtime;
results.pi0 =                   pi0;
results.pi0_res =               pi0_relRes;
results.pi0_flag =              pi0_flag;
results.pi0_time =              pi0_time;
results.unitCell_solvability =  unitCell_solvability;
results.unitCell_RHS =          unitCell_RHS;
results.unitCell_soln =         unitCell_soln;
results.unitCell_relRes =       unitCell_relRes;
results.unitCell_flag =         unitCell_flag;
results.unitCell_time =         unitCell_time;
results.Deff_mat =              Deff;
results.Deff_term1 =            DeffTerm1;
results.Deff_term2 =            DeffTerm2;
results.Deff =                  Deff(1);
results.geometry =              latticeGeo;
results.mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps );

%% slowOneSite
% geometry
dim = 2;
m = 2;
name = 'square';
obRad = 0;
obCtr = [.5 .5];
diagJumps = 0;
specialSetting = 'm2_slowOneSite'; %slowOneSite
driftMult = 0;
driftDecay = [];
obSlowdownFctr = .1;
bdyDist = [];

% Monte Carlo
numTraj = 100;
startNodeInd = [];
negateSteps = 1;
plotOn = 0;

latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);

pi0 = (1/4)*ones(4,1);
pi0_relRes = [];
pi0_flag = [];
pi0_time = [];

unitCell_soln = [pi0 pi0];
unitCell_RHS = [];
unitCell_relRes = [];
unitCell_solvability = [];
unitCell_flag = [];
unitCell_time = [];

LUtime = [];

[Deff, DeffTerm1, DeffTerm2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln );

results.L =                     L;
results.nodes =                 nodes;
results.edges =                 edges;
results.edgeRates =             edgeRates;
results.edgeJumps =             edgeJumps;
results.LU_time =               LUtime;
results.pi0 =                   pi0;
results.pi0_res =               pi0_relRes;
results.pi0_flag =              pi0_flag;
results.pi0_time =              pi0_time;
results.unitCell_solvability =  unitCell_solvability;
results.unitCell_RHS =          unitCell_RHS;
results.unitCell_soln =         unitCell_soln;
results.unitCell_relRes =       unitCell_relRes;
results.unitCell_flag =         unitCell_flag;
results.unitCell_time =         unitCell_time;
results.Deff_mat =              Deff;
results.Deff_term1 =            DeffTerm1;
results.Deff_term2 =            DeffTerm2;
results.Deff =                  Deff(1);
results.geometry =              latticeGeo;
results.mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps );