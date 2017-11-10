function [soln, RHS, relRes, solvability, flag, time] = unitCell(L,LTFactors,pi0,nodes,edges,edgeRates,edgeJumps,TOL,verbose)

if nargin < 8 || isempty(TOL)
    nNodes = size(L,1);
    TOL = 1e-8;sqrt(eps(1/nNodes));
end

if nargin < 9 || isempty(verbose)
    verbose = 1;
end

if verbose
    fprintf('\tUnit-cell problem:\n');
end
tic

LL = LTFactors.L; % unpacking LTFactors should not use memory
U = LTFactors.U;
P = LTFactors.P;
Q = LTFactors.Q;
R = LTFactors.R;

flag = 1;

% check solvability of unit-cell problem
solvability = sum(edgeJumps.*edgeRates.*pi0(edges(:,1)),1);
if any(solvability > 1e-10)
    warning('Unit cell problem may not be solvable. Solvability res = %.4e.',norm(solvability));
end

% calculate RHS of unit-cell problem
RHS = getUnitCellRHS(L,nodes,edges,edgeJumps, pi0);

% RHS can be zero. In this case, solution is pi0.
if norm(RHS(:)) < TOL
    dim = size(edgeJumps,2);
    soln = repmat(pi0,1,dim);
    
    err = L'*soln - RHS;
    normErr = sqrt(sum(err.^2,1));
    relRes = normErr;
else
    
    MSGID1 = 'MATLAB:nearlySingularMatrix';
    MSGID2 = 'MATLAB:singularMatrix';
    warning('off', MSGID1)
    warning('off', MSGID2)

    % solve via LU factors
    soln = Q * (U \ (LL \ (P * (R \ RHS))));

    err = L'*soln - RHS;
    normErr = sqrt(sum(err.^2,1));
    normRHS = sqrt(sum(RHS.^2,1));
    relRes = normErr./normRHS;

    % try backslash operator if LU backslash fails
    if any(relRes > TOL) || any(~isfinite(soln(:)))

        if verbose
            fprintf('\t\tLU Solver failed (rel res = %.4e). Trying backslash.\n', max(relRes));
        end
        soln = L'\RHS;

        err = L'*soln - RHS;
        normErr = sqrt(sum(err.^2,1));
        normRHS = sqrt(sum(RHS.^2,1));
        relRes = normErr./normRHS;

    end
    
    if any(~isfinite(soln(:)))
        if verbose
            fprintf('\t\tSolution has non-finite entries. Perturbing L.\n');
        end
        
        perturb = rand(nnz(L),1);
        inds = L ~= 0;
        Lperturb = L;
        Lperturb(inds) = Lperturb(inds) + perturb;
        soln = (Lperturb)'\RHS;

        err = L'*soln - RHS;
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
if verbose
    fprintf('\tCalculated unit-cell solution in %.1f seconds.\n',time);
end
    
end

function RHS = getUnitCellRHS(L,nodes,edges,edgeJumps, pi0) %L, pi0, nodes, edges, edgeJumps)
tic

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