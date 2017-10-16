function [L,nodes,edges,edgeRates,edgeJumps] = latticeSetup(ghParams)

m =         ghParams.m;
dim =       ghParams.dim;
R =         ghParams.R;
geometry =  ghParams.geometry;
ctr =       ghParams.ctr;

[nodes, nodeInds] = getNodes_lattice( R, m, dim, geometry, ctr );

[L,edges,edgeRates,edgeJumps] = getRateMat_lattice( nodes, nodeInds, ghParams );

if m == 2 && ~strcmpi(ghParams.specialSetting_m2,'none')
    
    edgeJumps(5:8,:) = -edgeJumps(5:8,:);
    edgeJumps(13:16,:) = -edgeJumps(13:16,:);

    if strcmpi(ghParams.specialSetting_m2,'blockOneSite')
        edgesToRemove = or(edges(:,1) == 4, edges(:,2) == 4);

        edges(edgesToRemove,:) = [];
        edgeJumps(edgesToRemove,:) = [];
        edgeRates(edgesToRemove) = [];
        nodes(4,:) = [];
        
    elseif strcmpi(ghParams.specialSetting_m2,'slowOneSite')
        edgesToSlow = or(edges(:,1) == 4, edges(:,2) == 4);
        delta = ghParams.rateCoeffs.delta;
        edgeRates(edgesToSlow) = delta*edgeRates(edgesToSlow);
    end
    L = sparse(edges(:,1),edges(:,2),edgeRates);
    L = L - diag(sum(L,2));
end

end