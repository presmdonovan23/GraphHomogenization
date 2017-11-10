function [factors, LUtime] = LUFull(L, verbose)

if nargin < 2 || isempty(verbose)
    verbose = 1;
end
tic;
[LL,U,P,Q,R] = lu(L);
LUtime = toc;

if verbose
    fprintf('\tCalculated LU decomposition in %.1f seconds.\n',LUtime);
end

factors.L = LL; % extra memory wont be allocated here
factors.U = U;
factors.P = P;
factors.Q = Q;
factors.R = R;

end