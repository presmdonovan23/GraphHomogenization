function [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo)

if ~isa(latticeGeo,'LatticeGeometry')
    error('Must supply a LatticeGeometry object.');
end
latticeGeo.validate;
if ~latticeGeo.isValid
    warning('Not a valid LatticeGeometry object. May have errors.');
end

[nodes, nodeInds] = getNodes_lattice( latticeGeo );

[L,edges,edgeRates,edgeJumps] = getRateMat_lattice(nodes,nodeInds,latticeGeo);

end

function [L,edges,edgeRates,edgeJumps] = getRateMat_lattice(nodes,nodeInds,latticeGeo)
% assumes lattice structure
%dim,geometryName,obRad,m,rateCoeffs,diagJumps,specialSetting_m2
tic

freeInds = nodeInds > 0;
nFree = sum(freeInds(:));

dim = latticeGeo.dim;
m = latticeGeo.m;
diagJumps = latticeGeo.diagJumps;

if dim == 2
    nbrShifts = [ 1  0; -1  0;...
                  0  1;  0 -1 ];
    if diagJumps
        nbrShifts = [nbrShifts;...
                     1  1;...
                     1 -1;...
                    -1  1;...
                    -1 -1 ];
    end
end

if dim == 3
    nbrShifts = [ 1  0  0;...
                 -1  0  0;...
                  0  1  0;...
                  0 -1  0;...
                  0  0  1;...
                  0  0 -1];
    if diagJumps
        nbrShifts = [nbrShifts;...
                     1  1  0;...
                     1 -1  0;...
                    -1  1  0;...
                    -1 -1  0;...
                     1  0  1;...
                     1  0 -1;...
                    -1  0  1;...
                    -1  0 -1;...
                     0  1  1;...
                     0  1 -1;...
                     0 -1  1;...
                     0 -1 -1];
    end
end

numNbrs = size(nbrShifts,1);

edgeStart = repmat(nodeInds(freeInds),numNbrs,1);
edgeEnd = zeros(nFree*numNbrs,1);
edgeRates = zeros(nFree*numNbrs,1);

for i = 1:numNbrs
    
    shift = nbrShifts(i,:);
    % shifting by [i,j,k] gives the neighbor: site + i*e_1 + j*e_2 + k*e_3
    nbrInds = circshift(nodeInds,shift);
    nonzeroRate = logical(nodeInds.*nbrInds);
    
    nonzeroRate = nonzeroRate(freeInds);
    
    nodeNbrs = nodes - latticeGeo.h*shift;
    curEdgeRates = rate_lattice(nodes,nodeNbrs,latticeGeo).*nonzeroRate;

    inds = 1 + (i-1)*nFree:i*nFree;
    edgeEnd(inds) = nbrInds(freeInds);
    edgeRates(inds) = curEdgeRates;
end

%% set up P

validInds = edgeRates > 0;
edgeStart = edgeStart(validInds);
edgeEnd = edgeEnd(validInds);
edgeRates = edgeRates(validInds);

L = sparse(edgeStart,edgeEnd,edgeRates);
L = L-diag(sum(L,2));

edges = [edgeStart,edgeEnd];

edgeJumps = nodes(edges(:,2),:) - nodes(edges(:,1),:);
inds = abs(edgeJumps) > .5;
edgeJumps(inds) = edgeJumps(inds) - sign(edgeJumps(inds));

%%
if m == 2 && ~strcmpi(specialSetting_m2,'none')
    
    edgeJumps(5:8,:) = -edgeJumps(5:8,:);
    edgeJumps(13:16,:) = -edgeJumps(13:16,:);

    if strcmpi(specialSetting_m2,'blockOneSite')
        edgesToRemove = or(edges(:,1) == 4, edges(:,2) == 4);

        edges(edgesToRemove,:) = [];
        edgeJumps(edgesToRemove,:) = [];
        edgeRates(edgesToRemove) = [];
        nodes(4,:) = [];
        
    elseif strcmpi(specialSetting_m2,'slowOneSite')
        edgesToSlow = or(edges(:,1) == 4, edges(:,2) == 4);
        delta = rateCoeffs.delta;
        edgeRates(edgesToSlow) = delta*edgeRates(edgesToSlow);
    end
    L = sparse(edges(:,1),edges(:,2),edgeRates);
    L = L - diag(sum(L,2));
end
%%
time = toc;
fprintf('Calculated rate matrix in %.1f seconds.\n',time);

end
