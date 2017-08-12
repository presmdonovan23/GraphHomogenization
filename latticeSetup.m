function [L,nodes,edges,edgeRates,edgeJumps,nodeInds] = latticeSetup(ghParams)

m =         ghParams.m;
dim =       ghParams.dim;
R =         ghParams.R;
geometry =  ghParams.geometry;
ctr =       ghParams.ctr;

[nodes, nodeInds] = getNodes_lattice( R, m, dim, geometry, ctr );

[L,edges,edgeRates,edgeJumps] = getRateMat_lattice( nodes, nodeInds, ghParams );

end