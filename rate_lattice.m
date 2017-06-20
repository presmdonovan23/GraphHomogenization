function val = rate_lattice(x,y,ghParams)
% Note that x = edge start and y = x + nu. This method will not work
% correctly if y - x != nu. E.g., [.1,.1] and [.1,.9].
% 
geometry = ghParams.geometry;
D0 = ghParams.D0;
h = ghParams.h;
nNodes = size(x,1);
alpha = ghParams.rateCoeffs.alpha;
diagJumps = ghParams.diagJumps;
K1 = ghParams.rateCoeffs.K1;
K2 = ghParams.rateCoeffs.K2;

if alpha ~= 0 && (diagJumps == 1 || strcmpi(geometry,'square'))
    error('Drift not yet implemented for diagonal jumps or square obstructions.')
end

if alpha == 0
    val = (D0/h^2)*ones(nNodes,1);
else % no diagonal jumps. all jumps along basis vectors
    ctr = .5;
    R = ghParams.R;
    
    jump = (abs(y - x) > 1e-12).*sign(y - x); % only 1 component will be nonzero
    mu = drift(x);
    
    d = sum(jump.*mu,2);  % contains the x,y, or z component of drift
    
    %val = D0/h^2 + alpha*d/(2*h);
    val = (1/h^2)*D0*exp(alpha*2*h*d);
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