function [Deff, term1, term2] = buildEffDiff( edges, edgeRates, edgeJumps, pi0, unitCell_soln, verbose )
% Computes effective diffusivity matrix given necessary ingredients
% Inputs:
%   1) edges = an e x 2 array where e = number of edges in the quotient
%   edge set. edges(i,1) = starting node index of the ith edge.
%   edges(i,2) = ending node index of the ith edge.
%   2) edgeRates = an e x 1 array where edgeRates(i) = jump rate of the
%   ith edge.
%   3) edgeJumps = an e x d array where edgeJumps(i) = the jump of the
%   ith edge. Equivalently, edgeJumps(i) = nodes(edges(i,2),:) -
%   nodes(edges(i,1),:).
%   4) pi0 = s x 1 array where pi0(i) = stationary distribution value at
%   ith node
%   5) unitCell_soln = s x 1 array where unitCell_soln(i) = unit-cell
%   solution at ith node
%   6) verbose = optional parameter controlling amount of output
% Outputs:
%   1) Deff = effective diffusivity matrix. Deff = term1 + term2.

if nargin < 6 || isempty(verbose)
    verbose = 1;
end

tic
dim = size(edgeJumps,2);

pi0_start = pi0(edges(:,1));
lambda_pi0 = edgeRates.*pi0_start;

omega_e = unitCell_soln(edges(:,1),:);

if dim == 2
    term1 = edgeJumps'*(edgeJumps.*[lambda_pi0 lambda_pi0]);
    term2 = edgeJumps'*(omega_e.*[edgeRates edgeRates]);
elseif dim == 3
    term1 = edgeJumps'*(edgeJumps.*[lambda_pi0 lambda_pi0 lambda_pi0]);
    term2 = edgeJumps'*(omega_e.*[edgeRates edgeRates edgeRates]);
end

term1 = .5*term1;
term2 = -.5*(term2 + term2');

Deff = term1 + term2;

time = toc;

if verbose
    fprintf('\tSet up effective diffusivity matrix in %.1f seconds.\n',time);
end

end