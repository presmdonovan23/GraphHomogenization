function [ghInput] = getDeff_whirlpool(ghInput,removeEdges,rateReduction,plotOn,plotLambda,plotNodeLabels)
% the purpose of this driver is to find a setting that satisfies unit-cell
% solvability but not detailed balance.
% we set up a simple setting with m = 3 and remove some edges to create a
% "cyclical" graph. this setting has some edges with 0 rates (or
% equivalently, e' is removed for some edges). we could also try making
% some rates very low if this creates issues.

%% modified ghInput
L = ghInput.L;

leftJumps = ismembertol(ghInput.edgeJumps,[-1/3,0],1e-10,'byrows',true);
rightJumps = ismembertol(ghInput.edgeJumps,[1/3,0],1e-10,'byrows',true);
downJumps = ismembertol(ghInput.edgeJumps,[0,-1/3],1e-10,'byrows',true);
upJumps = ismembertol(ghInput.edgeJumps,[0,1/3],1e-10,'byrows',true);

edgeStartNodes = ghInput.nodes(ghInput.edges(:,1),:);
edgeEndNodes = ghInput.nodes(ghInput.edges(:,2),:);

edgesToDamp1 = find(and(rightJumps,edgeStartNodes(:,2) > .5));   % edges on upper side that go right
edgesToDamp2 = find(and(leftJumps,edgeStartNodes(:,2) < .5));    % edges on lower side that go left
edgesToDamp3 = find(and(downJumps,edgeStartNodes(:,1) > .5));    % edges on right side that go down
edgesToDamp4 = find(and(upJumps,edgeStartNodes(:,1) < .5));      % edges on left side that go up

edgesToDamp = [edgesToDamp1; edgesToDamp2; edgesToDamp3; edgesToDamp4];
edgesToDamp = unique(edgesToDamp);

ghInput.edgeRates(edgesToDamp) = rateReduction*ghInput.edgeRates(edgesToDamp);
for i = 1:length(edgesToDamp)
    edge = ghInput.edges(edgesToDamp(i),:);
    L(edge(1),edge(2)) = rateReduction*L(edge(1),edge(2));
end

edgesToRemove1 = find(and(rightJumps,and(edgeStartNodes(:,1) > .6,edgeStartNodes(:,2) == .5)));  % edges on right center boundary that jump right
edgesToRemove2 = find(and(upJumps,and(edgeStartNodes(:,1) == .5,edgeStartNodes(:,2) > .6)));  % edges on upper center boundary that jump up
edgesToRemove3 = find(and(leftJumps,and(edgeStartNodes(:,1) < .4,edgeStartNodes(:,2) == .5)));  % edges on left center boundary that jump left
edgesToRemove4 = find(and(downJumps,and(edgeStartNodes(:,1) == .5,edgeStartNodes(:,2) < .4)));  % edges on lower center boundary that jump down

edgesToRemove = [edgesToRemove1; edgesToRemove2; edgesToRemove3; edgesToRemove4 ]; 
edgesToRemove = [edgesToRemove; edgesToDamp];
edgesToRemove = unique(edgesToRemove);
if ~removeEdges
    edgesToRemove = [];
end
for i = 1:length(edgesToRemove)
    edge = ghInput.edges(edgesToRemove(i),:);
    L(edge(1),edge(2)) = 0;
end
ghInput.edges(edgesToRemove,:) = [];
ghInput.edgeRates(edgesToRemove) = [];
ghInput.edgeJumps(edgesToRemove,:) = [];

L = L - diag(diag(L));
L = L - diag(sum(L,2));
s = 6/max(abs(L(:)));
L = L*s;

ghInput.L = L;
ghInput.edgeRates = ghInput.edgeRates*s;

end
