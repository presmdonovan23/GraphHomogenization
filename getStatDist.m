function [pi0, relRes, flag, time] = getStatDist(L,LTFactors, TOL)
% ** add some checks for how this might fail. maybe do several iterations
% and choose the best stationary distribution.
tic

nNodes = size(L,1);
if nargin < 3
    TOL = sqrt(eps(1/nNodes));
end

flag = 1;

LL = LTFactors.L; % unpacking LTFactors should not use memory
U = LTFactors.U;
P = LTFactors.P;
Q = LTFactors.Q;
R = LTFactors.R;

pi0 = randn(nNodes,1);
pi0 = pi0./max(abs(pi0));

MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID)
pi0 = Q * (U \ (LL \ (P * (R \ pi0))));
warning('on', MSGID)    

pi0 = pi0./sum(pi0);

relRes = norm(L'*pi0)/norm(pi0);

if relRes > TOL
    flag = 2;
    warning('Error exceeds tolerance (|L^T * pi0| = %.10f.',relRes)
end

time = toc;
fprintf('Calculated stationary distribution in %.1f seconds.\n',time);

end