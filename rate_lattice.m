function val = rate_lattice(x,y,ghParams)
% Note that x = edge start and y = x + nu. This method will not work
% correctly if y - x != nu. E.g., [.1,.1] and [.1,.9].
% You can assume the edge (x,y) exists and has a positive rate. This is
% taken care of in the parent function, getRateMat_lattice.
geometry = ghParams.geometry;
D0 = ghParams.D0;
h = ghParams.h;
nNodes = size(x,1);
alpha = ghParams.rateCoeffs.alpha;
diagJumps = ghParams.diagJumps;
K1 = ghParams.rateCoeffs.K1;
K2 = ghParams.rateCoeffs.K2;

if alpha ~= 0 && (diagJumps >= 1 || strcmpi(geometry,'square'))
    error('Drift not yet implemented for diagonal jumps or square obstructions.')
end

if diagJumps == 0
    lambda = D0/h^2;
else
    lambda = D0/(3*h^2);
end

if alpha == 0
    val = lambda*ones(nNodes,1);
    
    if ghParams.diagJumps == 2 && strcmpi(geometry,'square') && ghParams.dim == 2
        val = correctDiag(val,lambda,x,y,ghParams);
    elseif ghParams.diagJumps == 2
        error('Corrected diagonal jumps only implemented for 2d square with no drift.');
    end
else % no diagonal jumps. all jumps along basis vectors
    
    ctr = .5;
    R = ghParams.R;
    
    jump = (abs(y - x) > 1e-12).*sign(y - x); % only 1 component will be nonzero
    mu = drift(x);
    
    d = sum(jump.*mu,2);  % contains the x,y, or z component of drift
    
    %val = lambda + alpha*d/(2*h);
    val = lambda*exp(alpha*2*h*d);
end

function val = drift(x)
% only makes sense for circular obstuctions

relX = x - ctr;
relNorm = sqrt(sum(relX.^2,2));

dist2ob = relNorm - R;
xNormalized = relX./relNorm;

val = K1*xNormalized.*exp(-K2*dist2ob);
val(~isfinite(val)) = 0;

end

end

function val = correctDiag(val,lambda,x,y,ghParams)
    h = ghParams.h;
    TOL = 1e-10;
    nu = (1/h)*(y-x);
    S = ghParams.R*2;
    isBlocked = @(site) (abs(site(:,1) - .5) < S/2) & (abs(site(:,2) - .5) < S/2);

    % technically this returns true if obstructed.
    border =            (abs(x(:,1) - .5) < S/2+h) & (abs(x(:,2) - .5) < S/2+h);

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
    val(case1a & border) = 2*lambda;
    case1b = upperJump & ((upperRightBlocked & rightBlocked) | (upperLeftBlocked & leftBlocked));
    val(case1b & border) = 2*lambda;

    case2a = leftJump & ((upperLeftBlocked & upperBlocked) | (lowerLeftBlocked & lowerBlocked));
    val(case2a & border) = 2*lambda;
    case2b = lowerJump & ((lowerRightBlocked & rightBlocked) | (lowerLeftBlocked & leftBlocked));
    val(case2b & border) = 2*lambda;

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