function [Deff, term1, term2] = calcDeffMatrix( edges, edgeRates, edgeJumps, pi0, unitCell_soln )

tic
dim = size(edgeJumps,2);

pi0_start = pi0(edges(:,1));
lambda_pi0 = edgeRates.*pi0_start;

omega_e = unitCell_soln(edges(:,1),:);

if dim == 2
    term1 = edgeJumps'*(edgeJumps.*[lambda_pi0 lambda_pi0]);
    term2 = edgeJumps'*(omega_e.*[edgeRates edgeRates]);
elseif dim == 3
    term1 = edgeJumps'*(edgeJumps.*[lambda_pi0 lambda_pi0 lambda_pi0]);
    term2 = edgeJumps'*(omega_e.*[edgeRates edgeRates edgeRates]);
end

Deff = .5*(term1 - term2 - term2');

time = toc;
fprintf('Set up effective diffusivity matrix in %.1f seconds.\n',time);

end