function results = getDeff_homog(ghInput)

L = ghInput.L;
[LL,U,P,Q,R] = lu(L');
%[LL,U,P,Q,R] = lu(L' - eps*speye(nNodes));
LTFactors.L = LL; % extra memeory wont be allocated here
LTFactors.U = U;
LTFactors.P = P;
LTFactors.Q = Q;
LTFactors.R = R;
clear LL U P Q R;

[pi0, pi0_relRes, pi0_flag, pi0_time] = getStatDist(L,LTFactors);

unitCell_RHS = getUnitCellRHS(L, pi0, ghInput.nodes, ghInput.edges);

[unitCell_soln, unitCell_relRes, unitCell_flag, unitCell_time] = getUnitCellSoln(L,LTFactors,unitCell_RHS);

Deff = calcDeffMatrix( ghInput.edges, ghInput.edgeRates, ghInput.edgeJumps, pi0, unitCell_soln );

results.L = ghInput.L;
results.nodes = ghInput.nodes;
results.edges = ghInput.edges;
results.edgeRates = ghInput.edgeRates;
results.edgeJumps = ghInput.edgeJumps;
results.pi0 = pi0;
results.pi0_res = pi0_relRes;
results.pi0_flag = pi0_flag;
results.pi0_time = pi0_time;
results.unitCell_RHS = unitCell_RHS;
results.unitCell_soln = unitCell_soln;
results.unitCell_relRes = unitCell_relRes;
results.unitCell_flag = unitCell_flag;
results.unitCell_time = unitCell_time;
results.Deff_mat = Deff;
results.Deff = Deff(1);

end

