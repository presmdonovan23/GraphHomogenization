numTraj = 0;
specialSetting_m2 = 'blockOneSite';
specialSetting_m2 = 'slowOneSite';
delta = .1; % only required if setting = slowOneSite

% DO NOT CHANGE
rho = 0;
m = 2;
dim = 2;
geometry = 'square';
diagJumps = 0;
D0 = 1;
rate = [];
rateCoeffs.alpha = 0;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;
rateCoeffs.delta = delta;
ctr = [];
negateSteps = 1;

%%
params = LatticeParams(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps,ctr);

[L,nodes,edges,edgeRates,edgeJumps] = latticeSetup(params);

LUtime = [];

if strcmpi(specialSetting_m2,'blockOneSite')
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
elseif strcmpi(specialSetting_m2,'slowOneSite')
    
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
end

[Deff,DeffTerm1,DeffTerm2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln );

MC = getDeff_MC( L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps );

results.params = params;
results.L = L;
results.nodes = nodes;
results.edges = edges;
results.edgeRates = edgeRates;
results.edgeJumps = edgeJumps;
results.LU_time = LUtime;
results.pi0 = pi0;
results.pi0_res = pi0_relRes;
results.pi0_flag = pi0_flag;
results.pi0_time = pi0_time;
results.unitCell_solvability = unitCell_solvability;
results.unitCell_RHS = unitCell_RHS;
results.unitCell_soln = unitCell_soln;
results.unitCell_relRes = unitCell_relRes;
results.unitCell_flag = unitCell_flag;
results.unitCell_time = unitCell_time;
results.Deff_mat = Deff;
results.Deff_term1 = DeffTerm1;
results.Deff_term2 = DeffTerm2;
results.Deff = Deff(1);

results.monteCarlo = MC;