function [omega, relRes, flag, time] = getUnitCellSoln(L,LTFactors,RHS,TOL)
tic

LL = LTFactors.L; % unpacking LTFactors should not use memory
U = LTFactors.U;
P = LTFactors.P;
Q = LTFactors.Q;
R = LTFactors.R;

flag = 1;

if nargin < 4
    nNodes = size(L,1);
    TOL = sqrt(eps(1/nNodes));
end

MSGID = 'MATLAB:nearlySingularMatrix';
warning('off', MSGID)

omega = Q * (U \ (LL \ (P * (R \ RHS))));
%omega = L'\RHS;
%eta = zeros(size(this.RHS));
%for i = 1:this.dim
%    eta(:,i) = bicgstab(this.L',this.RHS(:,i),10^-8,2000);
%end
warning('on', MSGID)

err = L'*omega - RHS;
normErr = sqrt(sum(err.^2,1));
normRHS = sqrt(sum(RHS.^2,1));
relRes = normErr./normRHS;

if any(relRes > TOL)
    flag = 2;
    
    dim = size(RHS,2);
    if dim == 2
        warning('Relative residual exceeds tolerance (relres = [%.2e,%.2e]).',relRes);
    elseif dim == 3
        warning('Relative residual exceeds tolerance (relres = [%.2e,%.2e,%.2e]).',relRes);
    end
        
end

time = toc;
fprintf('Calculated unit-cell solution in %.1f seconds.\n',time);
    
end
