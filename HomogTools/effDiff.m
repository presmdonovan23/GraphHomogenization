function results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry,verbose)
% Inputs:
%   1) L = rate matrix
%   2) nodes = an s x d array where s = number of nodes in quotient node
%   set. s(i,1) = coordinates of ith node.
%   3) edges = an e x 2 array where e = number of edges in the quotient
%   edge set. edges(i,1) = starting node index of the ith edge.
%   edges(i,2) = ending node index of the ith edge.
%   4) edgeRates = an e x 1 array where edgeRates(i) = jump rate of the
%   ith edge.
%   5) edgeJumps = an e x d array where edgeJumps(i) = the jump of the
%   ith edge. Equivalently, edgeJumps(i) = nodes(edges(i,2),:) -
%   nodes(edges(i,1),:).
%   6) numTraj = number of Monte Carlo trajectories
%   7) startNodeInd = index of node that Monte Carlo simulations start at.
%   Can be empty, scalar (all trajectories start at same node), or numTraj
%   x 1 (one entry for the initial condition of each trajectory
%   8) geometry = latticeGeometry object or empty
%   9) verbose = optional parameter controlling amount of output
if nargin < 6 || isempty(numTraj)
    numTraj = 0;
end
if nargin < 7 || isempty(startNodeInd)
    if numTraj == 0
        startNodeInd = [];
    else
        startNodeInd = 1;
    end
end
if nargin < 8
    geometry = [];
end
if nargin < 9 || isempty(verbose)
    verbose = 1;
end

negateSteps = 0;

if isa(geometry,'LatticeGeometry')
    if geometry.m == 2
        negateSteps = 1;
    end
end

plotOn = 0;

% Calc Deff via homog theory
results = effDiff_homog(L,nodes,edges,edgeRates,edgeJumps,verbose);
% we story geometry object but do not need to use it at this point
results.geometry = geometry;
% Calc Deff via Monte Carlo
results.mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps, verbose );

end