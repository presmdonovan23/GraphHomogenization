function [ cond1, cond2, cond3 ] = checkDetailedBalance( nodes,edges,edgeRates,L, pi0 )
%CHECKDETAILEDBALANCE Summary of this function goes here
%   Detailed explanation goes here

%% check detailed balance
nNodes  = size(nodes,1);
nEdges = size(edges,1);
g = rand(nNodes,1);
f = rand(nNodes,1);

rates = edgeRates;

cond1 = 0;
cond2 = zeros(nEdges,1);

ratesPrime = zeros(nEdges,1);
for i = 1:nEdges
    % check L = L^T wrt pi
    edgeStart = edges(i,1);
    edgeEnd = edges(i,2);
    cond1 = cond1 + (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart);
    
    % check other thing
    ratesPrime(i) = full(L(edgeEnd,edgeStart));
    cond2(i) = (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart) + ...
        (g(edgeStart)*f(edgeEnd) - f(edgeStart)*g(edgeEnd))*ratesPrime(i)*pi0(edgeEnd);
    
end
cond3 = pi0(edges(:,1)).*edgeRates - pi0(edges(:,2)).*ratesPrime;

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

