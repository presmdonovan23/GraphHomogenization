function [L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(...
    dim,geometryName,obRad,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2)

if nargin < 7 || isempty(rateCoeffs)
    rateCoeffs.K1 = 0;
    rateCoeffs.K2 = 0;
    rateCoeffs.alpha = 0;
    rateCoeffs.delta = 0;
    rateCoeffs.dist = 0;
end
    
if nargin < 8 || isempty(diagJumps)
    diagJumps = 0;
end
if nargin < 9 || isempty(ctr)
    ctr = .5;
end
if nargin < 10 || isempty(specialSetting_m2)
    specialSetting_m2 = 'none';
end

[nodes, nodeInds] = getNodes_lattice( obRad, m, dim, geometryName, ctr );

[L,edges,edgeRates,edgeJumps] = getRateMat_lattice(nodes,nodeInds,dim,geometryName,obRad,m,rateCoeffs,diagJumps,ctr);

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

end

function [L,edges,edgeRates,edgeJumps] = getRateMat_lattice( ...
    nodes,nodeInds,dim,geometryName,obRad,m,rateCoeffs,diagJumps,ctr)
% assumes lattice structure

tic

h = 1/m;

freeInds = nodeInds > 0;
nFree = sum(freeInds(:));

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
    
    alpha = rateCoeffs.alpha;
    K1 = rateCoeffs.K1;
    K2 = rateCoeffs.K2;
    
    shift = nbrShifts(i,:);
    % shifting by [i,j,k] gives the neighbor: site + i*e_1 + j*e_2 + k*e_3
    nbrInds = circshift(nodeInds,shift);
    nonzeroRate = logical(nodeInds.*nbrInds);
    
    if diagJumps == 2
        isNbrFree{i} = nonzeroRate;
    end
    nonzeroRate = nonzeroRate(freeInds);
    
    nodeNbrs = nodes - h*shift;
    curEdgeRates = rate_lattice(nodes,nodeNbrs,dim,geometryName,obRad,h,alpha,diagJumps,K1,K2).*nonzeroRate;

    inds = 1 + (i-1)*nFree:i*nFree;
    edgeEnd(inds) = nbrInds(freeInds);
    edgeRates(inds) = curEdgeRates;
end

if strcmpi(geometryName,'circleSlowdown')
    % Slows all rates of edges starting or ending in square
    delta = rateCoeffs.delta;
    
    nodesInObs = find(sqrt(sum((nodes - ctr).^2,2)) <= R);
    slowedEdges = or(ismember(edgeStart,nodesInObs),ismember(edgeEnd,nodesInObs));
    edgeRates(slowedEdges) = edgeRates(slowedEdges)*delta;
elseif strcmpi(geometryName,'squareSlowdown')
    % Slows all rates of edges starting or ending in square
    delta = rateCoeffs.delta;
    
    relCoords = nodes - ctr;
    nodesInObs = find(all(abs(relCoords) < obRad,2));
    
    slowedEdges = or(ismember(edgeStart,nodesInObs),ismember(edgeEnd,nodesInObs));
    edgeRates(slowedEdges) = edgeRates(slowedEdges)*delta;
elseif strcmpi(geometryName,'squareBonding')
    % Slows all rates leaving nodes that are within dist
    delta = rateCoeffs.delta;
    dist = rateCoeffs.dist;
    
    relCoords = nodes - ctr;
    nodesAtBdy = find(all(abs(relCoords) < (obRad + dist),2));
    
    slowedEdges = ismember(edgeStart,nodesAtBdy);
    edgeRates(slowedEdges) = edgeRates(slowedEdges)*delta;
elseif strcmpi(geometryName,'squareBdyRepel')
    % Slows all rates from nodes at boundary to nodes not at boundary
    delta = rateCoeffs.delta;
    dist = rateCoeffs.dist;
    
    relCoords = nodes - ctr;
    nodesAtBdy = all(abs(relCoords) < (obRad + dist),2);
    nodesNotAtBdy = ~nodesAtBdy;
    
    nodesAtBdy = find(nodesAtBdy);
    nodesNotAtBdy = find(nodesNotAtBdy);
    
    acceleratedEdges = ismember(edgeStart,nodesAtBdy) & ismember(edgeEnd,nodesNotAtBdy);
    edgeRates(acceleratedEdges) = edgeRates(acceleratedEdges)*delta;
    
    slowedEdges = ismember(edgeStart,nodesNotAtBdy) & ismember(edgeEnd,nodesAtBdy);
    edgeRates(slowedEdges) = edgeRates(slowedEdges)/delta;
elseif strcmpi(geometryName,'squareBdyAttract')
    % Increase all rates from nodes not at boundary to nodes at boundary
    delta = rateCoeffs.delta;
    dist = rateCoeffs.dist;
    
    relCoords = nodes - ctr;
    nodesAtBdy = all(abs(relCoords) < (obRad + dist),2);
    nodesNotAtBdy = ~nodesAtBdy;
    
    nodesAtBdy = find(nodesAtBdy);
    nodesNotAtBdy = find(nodesNotAtBdy);
    
    acceleratedEdges = ismember(edgeStart,nodesNotAtBdy) & ismember(edgeEnd,nodesAtBdy);
    edgeRates(acceleratedEdges) = edgeRates(acceleratedEdges)*delta;
    
    slowedEdges = ismember(edgeStart,nodesAtBdy) & ismember(edgeEnd,nodesNotAtBdy);
    edgeRates(slowedEdges) = edgeRates(slowedEdges)/delta;
elseif strcmpi(geometryName,'squareBdySlow')
    % Slow all rates from boundary node to boundary node.
    % Meant to simulate nodes sitting on obstruction boundary.
    delta = rateCoeffs.delta;
    dist = rateCoeffs.dist;
    
    trueObRad = obRad + .25*h;
    relCoords = nodes - ctr;
    nodesAtBdy = all(abs(relCoords) < (obRad + dist),2);
    nodesNotAtBdy = ~nodesAtBdy;
    nodesAtCorner = all(abs(abs(relCoords) - trueObRad) < 1e-10,2);
    
    nodesAtBdy = find(nodesAtBdy);
    nodesNotAtBdy = find(nodesNotAtBdy);
    nodesAtCorner = find(nodesAtCorner);
    
    acceleratedEdges = ismember(edgeStart,nodesAtBdy) & ismember(edgeEnd,nodesNotAtBdy) & ~ismember(edgeStart,nodesAtCorner);
    edgeRates(acceleratedEdges) = edgeRates(acceleratedEdges)*delta;
    
    % need to do some surgery here because squareBdySlow is using a hack to
    % work.
    if diagJumps > 0
        warning('squareBdySlow geometry with diagonal jumps may not work as intended.');
        validInds = edgeRates > 0;
        edgeStart = edgeStart(validInds);
        edgeEnd = edgeEnd(validInds);
        edgeRates = edgeRates(validInds);
        edges = [edgeStart,edgeEnd];
        edgeJumps = nodes(edges(:,2),:) - nodes(edges(:,1),:);
        inds = abs(edgeJumps) > .5;
        edgeJumps(inds) = edgeJumps(inds) - sign(edgeJumps(inds));

        diagJumps = ismembertol(abs(edgeJumps),[h,h],1e-10,'byrows',1);
        bdyJumps = ismember(edges(:,1),nodesAtBdy) & ismember(edges(:,2),nodesAtBdy);
        edgeRates(diagJumps & bdyJumps) = 0;
        
        diagJumpsFromBdy = ismember(edgeStart,nodesAtBdy) & ismember(edgeEnd,nodesNotAtBdy) & diagJumps;
        edgeRates(diagJumpsFromBdy) = edgeRates(diagJumpsFromBdy)/delta;
    end
    %edges(diagJumps,:) = [];
    %edgeJumps(diagJumps,:) = [];
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

time = toc;
fprintf('Calculated rate matrix in %.1f seconds.\n',time);

end
