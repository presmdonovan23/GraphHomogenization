function results = getDeff_homog_D2_m2(ghInput, specialSetting_m2)

L = ghInput.L;
tic;
[LL,U,P,Q,R] = lu(L');
LUtime = toc;
fprintf('Calculated LU decomposition in %.1f seconds.\n',LUtime);

LTFactors.L = LL; % extra memory wont be allocated here
LTFactors.U = U;
LTFactors.P = P;
LTFactors.Q = Q;
LTFactors.R = R;
clear LL U P Q R;

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

[Deff,DeffTerm1,DeffTerm2] = calcDeffMatrix( ghInput.edges, ghInput.edgeRates, ghInput.edgeJumps, pi0, unitCell_soln );

results.L = ghInput.L;
results.nodes = ghInput.nodes;
results.edges = ghInput.edges;
results.edgeRates = ghInput.edgeRates;
results.edgeJumps = ghInput.edgeJumps;
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

end
