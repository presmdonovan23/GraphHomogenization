function results = effDiff_homog(L,nodes,edges,edgeRates,edgeJumps)

% get LU factors
[LTFactors, LUtime] = LUFull(L');
% solve for stationary distribution
[pi0, pi0_relRes, pi0_flag, pi0_time] = statDist(L,LTFactors);
% solve unit-cell problem
[unitCell_soln, unitCell_RHS, unitCell_relRes, unitCell_solvability, unitCell_flag, unitCell_time] = ...
    unitCell(L,LTFactors,pi0,nodes,edges,edgeRates,edgeJumps);
% calculate effective diffusivity
[Deff,DeffTerm1,DeffTerm2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln );

% store everything
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

end
