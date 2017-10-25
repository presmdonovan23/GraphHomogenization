function [factors, LUtime] = LUFull(L)
tic;
[LL,U,P,Q,R] = lu(L);
LUtime = toc;
fprintf('Calculated LU decomposition in %.1f seconds.\n',LUtime);

factors.L = LL; % extra memory wont be allocated here
factors.U = U;
factors.P = P;
factors.Q = Q;
factors.R = R;

end