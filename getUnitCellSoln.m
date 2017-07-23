function [omega, RHS, relRes, solvability, flag, time] = getUnitCellSoln(L,LTFactors,pi0,ghInput,TOL)
tic

LL = LTFactors.L; % unpacking LTFactors should not use memory
U = LTFactors.U;
P = LTFactors.P;
Q = LTFactors.Q;
R = LTFactors.R;

flag = 1;

if nargin < 5 || isempty(TOL)
    nNodes = size(L,1);
    TOL = sqrt(eps(1/nNodes));
end

solvability = sum(ghInput.edgeJumps.*ghInput.edgeRates.*pi0(ghInput.edges(:,1)),1);
if any(solvability > 1e-10)
    warning('Unit cell problem may not be solvable. Solvability res = %.4e.',norm(unitCell_solvability));
end

RHS = getUnitCellRHS(ghInput, pi0);

% RHS can be zero. In this case, solution is pi0.
if norm(RHS(:)) < TOL
    dim = size(ghInput.edgeJumps,2);
    omega = repmat(pi0,1,dim);
    
    err = L'*omega - RHS;
    normErr = sqrt(sum(err.^2,1));
    relRes = normErr;
else
    
    MSGID1 = 'MATLAB:nearlySingularMatrix';
    MSGID2 = 'MATLAB:singularMatrix';
    warning('off', MSGID1)
    warning('off', MSGID2)

    omega = Q * (U \ (LL \ (P * (R \ RHS))));

    err = L'*omega - RHS;
    normErr = sqrt(sum(err.^2,1));
    normRHS = sqrt(sum(RHS.^2,1));
    relRes = normErr./normRHS;

    if any(relRes > TOL) || any(~isfinite(omega(:)))

        warning('LU Solver failed. Trying backslash.');
        omega = L'\RHS;

        err = L'*omega - RHS;
        normErr = sqrt(sum(err.^2,1));
        normRHS = sqrt(sum(RHS.^2,1));
        relRes = normErr./normRHS;

    end

    if any(relRes > TOL) || any(~isfinite(omega(:)))

        warning('Backslash failed. Perturbing matrix.');
        omega = (L+eps*rand(size(L)))'\RHS;

        err = L'*omega - RHS;
        normErr = sqrt(sum(err.^2,1));
        normRHS = sqrt(sum(RHS.^2,1));
        relRes = normErr./normRHS;

    end

    warning('on', MSGID1)
    warning('on', MSGID2)
    
end

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

function RHS = getUnitCellRHS(ghInput, pi0) %L, pi0, nodes, edges, edgeJumps)
tic

L = ghInput.L;
nodes = ghInput.nodes;
edges = ghInput.edges;
edgeJumps = ghInput.edgeJumps;

[nNodes,dim] = size(nodes);

%inJumps = nodes(edges(:,1),:) - nodes(edges(:,2),:);
%inds = abs(inJumps) > .5;
%inJumps(inds) = inJumps(inds) - sign(inJumps(inds));
inJumps = -edgeJumps;

rows = edges(:,2);
cols = edges(:,1);
pi0_e = sparse(rows,cols,pi0(edges(:,1)));

RHS = zeros(nNodes,dim);
for i = 1:dim
    
    %rows = edges(:,1);
    %cols = edges(:,2);
    nu_ei = sparse(rows,cols,edgeJumps(:,i));
    
    RHS(:,i) = sum(pi0_e.*nu_ei.*L',2);

end

end