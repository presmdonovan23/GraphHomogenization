% the purpose of this driver is to find a setting that satisfies unit-cell
% solvability but not detailed balance.
% we set up a simple setting with m = 3 and remove some edges to create a
% "cyclical" graph. this setting has some edges with 0 rates (or
% equivalently, e' is removed for some edges). we could also try making
% some rates very low if this creates issues.
removeEdges = 0;
rateReduction = .5;
rho = .25;
m = 3;

dim = 2;
geometry = 'circle';
diagJumps = 0;
D0 = 1;
alpha = 0;
K1 = 25;
K2 = 10;


numTraj = 8000;
for i = 1:8
    startNodeInd(1 + (i-1)*1000:i*1000) = i;
end
%startNodeInd = 1;
plotOn = 0;
rate = []; % field is current obsolete
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;

h = 1/m;
        
ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

%% modified ghInput
L = ghInput.L;

leftJumps = ismembertol(ghInput.edgeJumps,[-1/m,0],1e-10,'byrows',true);
rightJumps = ismembertol(ghInput.edgeJumps,[1/m,0],1e-10,'byrows',true);
downJumps = ismembertol(ghInput.edgeJumps,[0,-1/m],1e-10,'byrows',true);
upJumps = ismembertol(ghInput.edgeJumps,[0,1/m],1e-10,'byrows',true);

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
    edgesToRemove = []; % delete this line to remove edges
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
%% get results
results_homog = getDeff_homog(ghInput)
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn )
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

term = 0;
ratesPrime = zeros(nEdges,1);
for i = 1:nEdges
    % check L = L^T wrt pi
    edgeStart = ghInput.edges(i,1);
    edgeEnd = ghInput.edges(i,2);
    term = term + (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart);
    
    % check other thing
    ratesPrime(i) = full(L(edgeEnd,edgeStart));
    term2(i) = (g(edgeEnd)*f(edgeStart) - f(edgeEnd)*g(edgeStart))*rates(i)*pi0(edgeStart) + ...
        (g(edgeStart)*f(edgeEnd) - f(edgeStart)*g(edgeEnd))*ratesPrime(i)*pi0(edgeEnd);
    
end

% check detailed balance
detailedBalance = results_homog.pi0(ghInput.edges(:,1)).*ghInput.edgeRates - results_homog.pi0(ghInput.edges(:,2)).*ratesPrime;


%% plot
%{
plotEdges = 1;
plotObs = 0;

nodes = ghInput.nodes;
edges = ghInput.edges;
edgeRates = ghInput.edgeRates;
edgeJumps = ghInput.edgeJumps;

dim = ghParams.dim;
geometry = ghParams.geometry;
obRad = ghParams.R;
maxHeadSize = .35*max(sqrt(sum(edgeJumps.^2,2)));

fh = figure;
hold on

plot(nodes(:,1),nodes(:,2),'b.','markersize',50);

if plotEdges
    quiver( nodes(edges(:,1),1),nodes(edges(:,1),2),...
            edgeJumps(:,1),edgeJumps(:,2),...
            'color','black',...
            'autoscale','off',...
            'MaxHeadSize',maxHeadSize,...
            'linewidth',3);
end
if plotObs

    cornerX = .5 - obRad;
    cornerY = .5 - obRad;
    pos = [cornerX cornerY 2*obRad 2*obRad];

    if strcmpi(geometry,'circle')
        curvature = [1 1];
    elseif strcmpi(geometry,'square')
        curvature = [0 0];
    end

    rectangle(  'Position', pos,...
                'facecolor','red',...
                'Curvature',curvature);

end


for i = 1:size(edges,1)
    edge = edges(i,:);
    rate = edgeRates(i);
    jump = edgeJumps(i,:);

    % plot rates as text
    textloc = nodes(edge(1),:) + .5*jump;
    %rateStr = sprintf('%d',round(rate));%num2str(rate,'%.2f')
    rateStr = '$\bar{\lambda}$';
    if jump(1) > 10e-6 %edge goes left to right
        textloc(2) = textloc(2) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');

    elseif jump(1) < -10e-6 %edge goes right to left
        textloc(2) = textloc(2) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');
        
        %set(h, 'rotation', 180)
    elseif jump(2) > 10e-6 %edge goes bottom to top
        textloc(1) = textloc(1) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');
        
        %set(h, 'rotation', 90)
    elseif jump(2) < -10e-6 %edge goes top to bottom
        textloc(1) = textloc(1) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');
        
        %set(h, 'rotation', 270)
    end
end

lh = legend('Nodes','Edges','location','southeast');
set(lh,'fontsize',18);
axis square
axis off
rectangle('position',[0 0 1 1],'linewidth',1)
%mySaveFig('detailedBalance',fh);
%}