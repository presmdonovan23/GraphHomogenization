function [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo, verbose)

if nargin < 2 || isempty(verbose)
    verbose = 1;
end

if verbose
    fprintf('Input setup (lattice):\n')
end
if ~isa(latticeGeo,'LatticeGeometry')
    error('Must supply a LatticeGeometry object.');
end

% check if LatticeGeometry object is valid
latticeGeo.validate;
if ~latticeGeo.isValid
    warning('Not a valid LatticeGeometry object. May have errors.');
end

[nodes, nodeInds] = getNodes_lattice( latticeGeo, verbose );

[L,edges,edgeRates,edgeJumps] = getRateMat_lattice(nodes,nodeInds,latticeGeo,verbose);

end

function [L,edges,edgeRates,edgeJumps] = getRateMat_lattice(nodes,nodeInds,latticeGeo,verbose)
% assumes lattice structure
%dim,geometryName,obRad,m,rateCoeffs,diagJumps,specialSetting_m2

if nargin < 4 || isempty(verbose)
    verbose = 1;
end

tic

freeInds = nodeInds > 0;
nFree = sum(freeInds(:));

dim = latticeGeo.dim;
m = latticeGeo.m;
diagJumps = latticeGeo.diagJumps;

% get possible jumps in 2D
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
% get possible jumps in 3D
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
    % in each iteration of this loop, we calculate all of the jump rates
    % of edges whose jump size is proportional nbrShfits(i,:)
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

%% set up L

validInds = edgeRates > 0;
edgeStart = edgeStart(validInds);
edgeEnd = edgeEnd(validInds);
edgeRates = edgeRates(validInds);

L = sparse(edgeStart,edgeEnd,edgeRates);
L = L - diag(sum(L,2));

edges = [edgeStart,edgeEnd];

% any edges with jump size > .5 must be an edge leaving the periodic cell
edgeJumps = nodes(edges(:,2),:) - nodes(edges(:,1),:);
inds = abs(edgeJumps) > .5;
edgeJumps(inds) = edgeJumps(inds) - sign(edgeJumps(inds));

%%
% if m = 2, we use a hard-coded fix because some edges will be incorrect.
% This is due to the fact that the code makes an assumption about the
% edges. See README.
specialSetting = latticeGeo.specialSetting;
if m == 2
    
    if size(edges,1) == 8
        edgesToNegate = [3,4,7,8];
    elseif size(edges,1) == 16
        edgesToNegate = [5:8 13:16];
    else
        error('Something went wrong. There should be 12 or 16 edges when m = 2.');
    end
    
    edgeJumps(edgesToNegate,:) = -edgeJumps(edgesToNegate,:);
    
    if strcmpi(specialSetting,'m2_slowOneSite')
        edgesToSlow = or(edges(:,1) == 4, edges(:,2) == 4);
        edgeRates(edgesToSlow) = latticeGeo.obSlowdownFctr*edgeRates(edgesToSlow);
    end
    
    L = sparse(edges(:,1),edges(:,2),edgeRates);
    L = L - diag(sum(L,2));
end
%%
time = toc;
if verbose
    fprintf('\tCalculated rate matrix in %.1f seconds.\n',time);
end

end
