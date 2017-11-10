function results = effDiff_homog(L,nodes,edges,edgeRates,edgeJumps,verbose)

if nargin < 6 || isempty(verbose)
    verbose = 1;
end
if verbose
    fprintf('Homogenization theory:\n');
end
TOL = [];

% get LU factors
[LTFactors, LUtime] = LUFull(L',verbose);
% solve for stationary distribution
[pi0, pi0_relRes, pi0_flag, pi0_time] = statDist(L,LTFactors,TOL,verbose);
% solve unit-cell problem
[unitCell_soln, unitCell_RHS, unitCell_relRes, unitCell_solvability, unitCell_flag, unitCell_time] = ...
    unitCell(L,LTFactors,pi0,nodes,edges,edgeRates,edgeJumps,TOL,verbose);
% calculate effective diffusivity (homog + monte carlo)
[Deff,DeffTerm1,DeffTerm2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln, verbose );

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
