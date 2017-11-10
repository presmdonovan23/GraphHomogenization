function [pi0, relRes, flag, time] = statDist(L, LTFactors, TOL, verbose)
% ** add some checks for how this might fail. maybe do several iterations
% and choose the best stationary distribution.
tic

nNodes = size(L,1);
if nargin < 3 || isempty(TOL)
    TOL = 1e-8;%sqrt(eps(1/nNodes));
end
if nargin < 4 || isempty(verbose)
    verbose = 1;
end

if verbose
    fprintf('\tStationary distribution:\n')
end
flag = 1;

LL = LTFactors.L; % unpacking LTFactors should not use memory
U = LTFactors.U;
P = LTFactors.P;
Q = LTFactors.Q;
R = LTFactors.R;

% initial guess for power method
pi0 = randn(nNodes,1);
pi0 = pi0./max(abs(pi0));

MSGID1 = 'MATLAB:nearlySingularMatrix';
MSGID2 = 'MATLAB:singularMatrix';
warning('off', MSGID1)
warning('off', MSGID2)

% calculate pi0 via LU factors (one iteration of power method)
pi0 = Q * (U \ (LL \ (P * (R \ pi0))));
pi0 = pi0./sum(pi0);

relRes = norm(L'*pi0)/norm(pi0);

% if error is bad, try another method
if relRes > TOL || any(~isfinite(pi0))
    
    if verbose
        fprintf('\t\tPower method failed (rel res = %.4e). Trying eigs.\n',relRes);
    end
    
    try
        % call MATLAB's eigs
        [pi0,~] = eigs(L',1,0);
    catch ME
        if verbose
            fprintf('\t\teigs failed. Trying null(full(L)) because eigs failed: %s\n',ME.message);
        end
        % last resort: convert matrix to full matrix and call null
        pi0 = null(full(L'));
    end
    
    pi0 = pi0./sum(pi0);
    relRes = norm(L'*pi0)/norm(pi0);

end

warning('on', MSGID1)
warning('on', MSGID2)

if relRes > TOL
    flag = 2;
    warning('Error exceeds tolerance (rel res = %.4e).',relRes)
end

time = toc;
if verbose
    fprintf('\t\tCalculated stationary distribution in %.1f seconds.\n',time);
end

end