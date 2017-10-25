removeEdges = 1;
rateReduction = .5;
numTraj = 00;
plotOn = 1;
plotLambda = 0;
plotNodeLabels = 1;
% DO NOT CHANGE
rho = .25;
m = 3;
dim = 2;
geometry = 'circle';
diagJumps = 0;
D0 = 1;
alpha = 0;
K1 = 25;
K2 = 10;
rate = [];
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;

for i = 1:8
    startNodeInd(1 + (i-1)*numTraj:i*numTraj) = i;
end

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(...
    dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);

%% need to make some modifications
leftJumps = ismembertol(edgeJumps,[-1/3,0],1e-10,'byrows',true);
rightJumps = ismembertol(edgeJumps,[1/3,0],1e-10,'byrows',true);
downJumps = ismembertol(edgeJumps,[0,-1/3],1e-10,'byrows',true);
upJumps = ismembertol(edgeJumps,[0,1/3],1e-10,'byrows',true);

edgeStartNodes = nodes(edges(:,1),:);
edgeEndNodes = nodes(edges(:,2),:);

edgesToDamp1 = find(and(rightJumps,edgeStartNodes(:,2) > .5));   % edges on upper side that go right
edgesToDamp2 = find(and(leftJumps,edgeStartNodes(:,2) < .5));    % edges on lower side that go left
edgesToDamp3 = find(and(downJumps,edgeStartNodes(:,1) > .5));    % edges on right side that go down
edgesToDamp4 = find(and(upJumps,edgeStartNodes(:,1) < .5));      % edges on left side that go up

edgesToDamp = [edgesToDamp1; edgesToDamp2; edgesToDamp3; edgesToDamp4];
edgesToDamp = unique(edgesToDamp);

edgeRates(edgesToDamp) = rateReduction*edgeRates(edgesToDamp);
for i = 1:length(edgesToDamp)
    edge = edges(edgesToDamp(i),:);
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
    edge = edges(edgesToRemove(i),:);
    L(edge(1),edge(2)) = 0;
end
edges(edgesToRemove,:) = [];
edgeRates(edgesToRemove) = [];
edgeJumps(edgesToRemove,:) = [];

L = L - diag(diag(L));
L = L - diag(sum(L,2));
s = 6/max(abs(L(:)));
L = L*s;

edgeRates = edgeRates*s;
%%

results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,geometry);

%%

if plotOn
    
    plotEdges = 1;
    plotObs = 0;
    
    obRad = rho/2;
    fh = drawWhirlpoolSetting( nodes,edges,edgeRates,edgeJumps,geometry,obRad, plotEdges, plotObs, plotLambda );

    % dirname = '';
    %mySaveFig([dirname 'detailedBalance2'],fh,'fig')
    %mySaveFig([dirname 'detailedBalance2'],fh,'png')
    %mySaveFig([dirname 'detailedBalance2'],fh,'eps')
    
end
