function results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry)

if nargin < 6 || isempty(numTraj)
    numTraj = 0;
end
if nargin < 8
    geometry = [];
end
plotOn = 0;

results = effDiff_homog(L,nodes,edges,edgeRates,edgeJumps);
results.geometry = geometry;
results.mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn );

end