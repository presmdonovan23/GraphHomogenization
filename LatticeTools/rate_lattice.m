function [rate,mu] = rate_lattice(startNode,endNode,latticeGeo)
% Note that x = edge start and y = x + nu. This method will not work
% correctly if y - x != nu. E.g., [.1,.1] and [.1,.9].
% You can assume the edge (x,y) exists and has a positive rate. This is
% taken care of in the parent function, getRateMat_lattice.

h = latticeGeo.h;
name = latticeGeo.name;
diagJumps = latticeGeo.diagJumps;
obCtr = latticeGeo.obCtr;
obRad = latticeGeo.obRad;
driftMult = latticeGeo.driftMult;
driftDecay = latticeGeo.driftDecay;
specialSetting = latticeGeo.specialSetting;
bdyDist = latticeGeo.bdyDist;

if diagJumps == 0
    lambda = 1/h^2;
else
    lambda = 1/(3*h^2);
end

if driftMult == 0
    nNodes = size(startNode,1);
    rate = lambda*ones(nNodes,1);
    
    if diagJumps == 2
        rate = correctDiag(rate,lambda,startNode,endNode,h,latticeGeo.sideLen);
    end
else % no diagonal jumps. all jumps along basis vectors
    
    jump = (abs(endNode - startNode) > 1e-12).*sign(endNode - startNode); % only 1 component will be nonzero
    mu = drift(startNode,obCtr,obRad,driftMult,driftDecay);
        
    d = sum(jump.*mu,2);  % contains the x,y, or z component of drift
    
    rate = lambda + d/(2*h);
    %rate = lambda*exp(driftMult*2*h*d);

end

if ~strcmpi(specialSetting,'none') && ~isempty(bdyDist)
    obSlowdownFctr = latticeGeo.obSlowdownFctr;
    bdyDist = latticeGeo.bdyDist;
    startNodeRel = startNode - obCtr;
    endNodeRel = endNode - obCtr;
    
    if strcmpi(name,'circle')    
        
        edgeStartInObs = sqrt(sum(startNodeRel.^2,2)) <= obRad;
        edgeEndInObs = sqrt(sum(endNodeRel.^2,2)) <= obRad;
        
    elseif strcmpi(name,'square')

        edgeStartInObs = all(abs(startNodeRel) < obRad,2);
        edgeEndInObs = all(abs(endNodeRel) < obRad,2);
        
        edgeStartInBdy = all(abs(startNodeRel) < (obRad + bdyDist),2);
        edgeEndInBdy = all(abs(endNodeRel) < (obRad + bdyDist),2);

    end
    
end

if strcmpi(specialSetting,'slowdown')
    % Slows all rates of edges starting or ending in square
    
    acceleratedEdges = [];
    slowedEdges = edgeStartInObs | edgeEndInObs;
    
elseif strcmpi(specialSetting,'bdyBonding')
    % Slows all rates leaving nodes that are within dis
    
    acceleratedEdges = [];
    slowedEdges = edgeStartInBdy;
    
elseif strcmpi(specialSetting,'bdyRepel')
    % Slows all rates from nodes at boundary to nodes not at boundary
    
    acceleratedEdges = edgeStartInBdy | ~edgeEndInBdy;
    slowedEdges = ~edgeStartInBdy | edgeEndInBdy;
    
elseif strcmpi(specialSetting,'bdyAttract')
    % Increase all rates from nodes not at boundary to nodes at boundary
    
    acceleratedEdges = ~edgeStartInBdy | edgeEndInBdy;
    slowedEdges = edgeStartInBdy | ~edgeEndInBdy;
    
elseif strcmpi(specialSetting,'bdySlow')
    % Slow all rates from boundary node to boundary node.
    % Meant to simulate nodes sitting on obstruction boundary.
    
    trueObRad = obRad + .25*h;
    
    edgeStartInBdy = all(abs(startNodeRel) < (trueObRad + bdyDist),2); % old code used obRad, not trueObRad
    edgeEndInBdy = all(abs(endNodeRel) < (trueObRad + bdyDist),2); % old code used obRad, not trueObRad
    edgeStartAtCorner = all(abs(abs(startNodeRel) - trueObRad) < 1e-10,2);
    
    acceleratedEdges = [];
    slowedEdges = edgeStartInBdy & ~edgeEndInBdy & ~edgeStartAtCorner;
    
    % need to do some surgery here because squareBdySlow is using a hack to work.
    if diagJumps > 0
        warning('squareBdySlow geometry with diagonal jumps may not work as intended.');

        edgeJumps = endNode - startNode;
        inds = abs(edgeJumps) > .5;
        edgeJumps(inds) = edgeJumps(inds) - sign(edgeJumps(inds));

        diagJumps = ismembertol(abs(edgeJumps),[h,h],1e-10,'byrows',1);
        bdyJumps = edgeStartInBdy & edgeEndInBdy;
        rate(diagJumps & bdyJumps) = 0;

        acceleratedEdges = edgeStartInBdy & ~edgeEndInBdy & diagJumps;
        
    end
end

if ~strcmpi(specialSetting,'none') && ~isempty(bdyDist)
    rate(slowedEdges) = rate(slowedEdges)*obSlowdownFctr;
    rate(acceleratedEdges) = rate(acceleratedEdges)/obSlowdownFctr;
end

end

function val = drift(x,obCtr,obRad,driftMult,driftDecay)
% only makes sense for circular obstuctions

relX = x - obCtr;
relNorm = sqrt(sum(relX.^2,2));

dist2ob = relNorm - obRad;
xNormalized = relX./relNorm;

val = driftMult*xNormalized.*exp(-driftDecay*dist2ob);
val(~isfinite(val)) = 0;

end

function val = correctDiag(val,lambda,x,y,h,sideLen)
    
    TOL = 1e-10;
    nu = (1/h)*(y-x);
    
    isBlocked = @(site) (abs(site(:,1) - .5) < sideLen/2) & (abs(site(:,2) - .5) < sideLen/2);

    % technically this returns true if obstructed.
    border =            (abs(x(:,1) - .5) < sideLen/2+h) & (abs(x(:,2) - .5) < sideLen/2+h);

    rightJump =         sum(abs(nu - [1, 0]),2) < TOL;
    leftJump =          sum(abs(nu - [-1, 0]),2) < TOL;
    upperJump =         sum(abs(nu - [0, 1]),2) < TOL;
    lowerJump =         sum(abs(nu - [0, -1]),2) < TOL;
    upperRightJump =    sum(abs(nu - [1, 1]),2) < TOL;
    upperLeftJump =     sum(abs(nu - [-1, 1]),2) < TOL;
    lowerRightJump =    sum(abs(nu - [1, -1]),2) < TOL;
    lowerLeftJump =     sum(abs(nu - [-1, -1]),2) < TOL;

    leftBlocked = isBlocked(x + h*[-1 0]);
    rightBlocked = isBlocked(x + h*[1 0]);
    lowerBlocked = isBlocked(x + h*[0 -1]);
    upperBlocked = isBlocked(x + h*[0 1]);
    upperRightBlocked = isBlocked(x + h*[1 1]);
    upperLeftBlocked =  isBlocked(x + h*[-1 1]);
    lowerRightBlocked = isBlocked(x + h*[1 -1]);
    lowerLeftBlocked =  isBlocked(x + h*[-1 -1]);

    case1a = rightJump & ((upperRightBlocked & upperBlocked) | (lowerRightBlocked & lowerBlocked));
    val(case1a & border) = 2*lambda; % change this to 1.5*lambda for a different but reasonable approach
    case1b = upperJump & ((upperRightBlocked & rightBlocked) | (upperLeftBlocked & leftBlocked));
    val(case1b & border) = 2*lambda; % change this to 1.5*lambda for a different but reasonable approach

    case2a = leftJump & ((upperLeftBlocked & upperBlocked) | (lowerLeftBlocked & lowerBlocked));
    val(case2a & border) = 2*lambda; % change this to 1.5*lambda for a different but reasonable approach
    case2b = lowerJump & ((lowerRightBlocked & rightBlocked) | (lowerLeftBlocked & leftBlocked));
    val(case2b & border) = 2*lambda; % change this to 1.5*lambda for a different but reasonable approach

    case3a = (upperJump | rightJump) & upperRightBlocked & ~upperBlocked & ~rightBlocked ;
    val(case3a & border) = 1.5*lambda;
    case3b = (upperJump | leftJump) & upperLeftBlocked & ~upperBlocked & ~leftBlocked ;
    val(case3b & border) = 1.5*lambda;
    case3c = (lowerJump | rightJump) & lowerRightBlocked & ~lowerBlocked & ~rightBlocked ;
    val(case3c & border) = 1.5*lambda;
    case3d = (lowerJump | leftJump) & lowerLeftBlocked & ~lowerBlocked & ~leftBlocked ;
    val(case3d & border) = 1.5*lambda;

    case4a = upperRightJump & ((upperBlocked & ~rightBlocked) | (~upperBlocked & rightBlocked));
    val(case4a & border) = .5*lambda;
    case4b = upperLeftJump & ((upperBlocked & ~leftBlocked) | (~upperBlocked  & leftBlocked));
    val(case4b & border) = .5*lambda;
    case4c = lowerRightJump & ((lowerBlocked & ~rightBlocked) | (~lowerBlocked  & rightBlocked));
    val(case4c & border) = .5*lambda;
    case4d = lowerLeftJump & ((lowerBlocked & ~leftBlocked) | (~lowerBlocked  & leftBlocked));
    val(case4d & border) = .5*lambda;
    
    case5a = upperJump & ((~upperRightBlocked & rightBlocked) | (~upperLeftBlocked & leftBlocked));
    val(case5a & border) = 1.5*lambda;
    case5b = lowerJump & ((~lowerRightBlocked & rightBlocked) | (~lowerLeftBlocked & leftBlocked));
    val(case5b & border) = 1.5*lambda;
    case5c = rightJump & ((~lowerRightBlocked & lowerBlocked) | (~upperRightBlocked & upperBlocked));
    val(case5c & border) = 1.5*lambda;
    case5d = leftJump & ((~lowerLeftBlocked & lowerBlocked) | (~upperLeftBlocked & upperBlocked));
    val(case5d & border) = 1.5*lambda;

end