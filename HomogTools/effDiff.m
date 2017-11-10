function results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry)

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

if isa(geometry,'LatticeGeometry')
    if geometry.m == 2
        negateSteps = 1;
    end
end

plotOn = 0;

results = effDiff_homog(L,nodes,edges,edgeRates,edgeJumps);
results.geometry = geometry;
results.mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps );

end