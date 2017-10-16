function [ cond1, cond2 ] = checkDetailedBalance( ghInput, results_homog )
%CHECKDETAILEDBALANCE Summary of this function goes here
%   Detailed explanation goes here

%% check detailed balance
nNodes  = size(ghInput.nodes,1);
nEdges = size(ghInput.edges,1);
g = rand(nNodes,1);
f = rand(nNodes,1);

rates = ghInput.edgeRates;
L = ghInput.L;
pi0 = results_homog.pi0;
pi0Start = pi0(ghInput.edges(:,1));
pi0End = pi0(ghInput.edges(:,2));

cond1 = 0;
ratesPrime = zeros(nEdges,1);
for i = 1:nEdges
    % check L = L^T wrt pi
    edgeStart = ghInput.edges(i,1);
    edgeEnd = ghInput.edges(i,2);
    cond1 = cond1 + (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart);
    
    % check other thing
    ratesPrime(i) = full(L(edgeEnd,edgeStart));
    cond2(i) = (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart) + ...
        (g(edgeStart)*f(edgeEnd) - f(edgeStart)*g(edgeEnd))*ratesPrime(i)*pi0(edgeEnd);
    
end
cond3 = results_homog.pi0(ghInput.edges(:,1)).*ghInput.edgeRates - results_homog.pi0(ghInput.edges(:,2)).*ratesPrime;

fprintf('\n');
fprintf('Unit cell solvability: %.10f.\n',sum(abs(results_homog.unitCell_solvability)));

if sum(abs(cond1)) < 1e-10
    fprintf('L = L^T wrt pi is satisfied (%.10f).\n',sum(abs(cond1)));
else
    fprintf('L = L^T wrt pi is not satisfied (%.10f).\n',sum(abs(cond1)));
end
% check detailed balance

if sum(abs(cond3)) < 1e-10
    fprintf('Detailed balance is satisfied (%.10f).\n',sum(abs(cond3)));
else
    fprintf('Detailed balance is not satisfied (%.10f).\n',sum(abs(cond3)));
end

end

