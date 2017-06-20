function RHS = getUnitCellRHS(L, pi0, nodes, edges)
tic

[nNodes,dim] = size(nodes);

inJumps = nodes(edges(:,1),:) - nodes(edges(:,2),:);
inds = abs(inJumps) > .5;
inJumps(inds) = inJumps(inds) - sign(inJumps(inds));

rows = edges(:,1);
cols = edges(:,2);
pi0_e = sparse(rows,cols,pi0(edges(:,2)));

RHS = zeros(nNodes,dim);
for i = 1:dim
    
    nu_ei = sparse(rows,cols,inJumps(:,i));
    
    RHS(:,i) = sum(pi0_e.*nu_ei.*L',2);

end

end