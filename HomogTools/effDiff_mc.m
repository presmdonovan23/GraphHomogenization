function mc = effDiff_mc(L, nodes, edges, edgeJumps, numTraj, startNodeInd, plotOn, negateSteps )

if nargin < 5 || isempty(numTraj)
    numTraj = 1000;
end
if nargin < 6 || isempty(startNodeInd)
    startNodeInd = 1;
end
if nargin < 7 || isempty(plotOn)
    plotOn = 0;
end
if nargin < 8 || isempty(negateSteps)
    negateSteps = 0;
end

if size(nodes,1) > 4 && negateSteps
    warning('negateSteps should only be 1 if m = 2.');
elseif size(nodes,1) <= 4 && ~negateSteps
    warning('It appears that m = 2 and negateSteps = 0. Results may be inaccurate if m = 2.')
end
    
if numTraj == 0
    mc.Deff = [];
    mc.Deff_var = [];
    mc.Deff_95CI = [];
    mc.Deff_all = [];
    mc.Deff_var_all = [];
    mc.numTraj = numTraj;
    mc.timeRec = [];
    mc.trajectoryPos = [];
    mc.startNodeInd = [];

    return;
end

tic

if mod(numTraj,100) ~= 0 && numTraj > 0    
    numTraj = 100*ceil(numTraj/100);
    fprintf('Rounding numTraj to %d.\n',numTraj);
end

if length(startNodeInd) == 1
    startNodeInd = repmat(startNodeInd,numTraj,1);
end

dim = size(edgeJumps,2);
numTimeRec = 50;
tmax = 25;

timeRec = linspace(0,tmax,numTimeRec);
locInterp = zeros(numTraj,dim,numTimeRec);


for p = 1:100
    for i = (1+(p-1)*numTraj/100):(p*numTraj/100)
        
        curStartNodeInd = startNodeInd(i);
        curStartNode = nodes(curStartNodeInd,:);
        
        [ timeCur, locCur ] = ...
            trajectory( L, edges, edgeJumps, curStartNode, curStartNodeInd, tmax, negateSteps );
        locInterp(i,:,:) = interp1(timeCur',locCur',timeRec')';
        
    end
    
    fprintf('%d%%,',p)
    if mod(p,10) == 0
        fprintf('\n');
    end
    
end

fprintf('\n');

sdInterp = squeeze(sum((locInterp - nodes(startNodeInd,:)).^2,2));
MSD = mean(sdInterp,1);

Deff_all = MSD(2:end)./(2*dim*timeRec(2:end)); % for all Deff estimates
Deff = Deff_all(end);

DeffVar_all = var(sdInterp./(2*dim*timeRec),[],1);
DeffVar = DeffVar_all(end);%var(sdInterp(:,end)/(2*dim*tmax));

Deff_CI_low = Deff - 1.96*sqrt(DeffVar)/sqrt(numTraj);
Deff_CI_high = Deff + 1.96*sqrt(DeffVar)/sqrt(numTraj);

time = toc;

mc.numTraj = numTraj;
mc.timeRec = timeRec;
mc.trajectoryPos = locInterp;
mc.startNodeInd = startNodeInd;
mc.Deff = Deff;
mc.Deff_all = Deff_all;
mc.Deff_var = DeffVar;
mc.Deff_var_all = DeffVar_all;
mc.Deff_95CI = [Deff_CI_low,Deff_CI_high];

if plotOn
    figure
    hold on
    plot(timeRec,sdInterp(1:15,:))
    plot(timeRec,MSD,'k-','linewidth',5)
    
    title('15 sample trajectories');
    xlabel('Time')
    ylabel('SD')
end

fprintf('Finished %d Monte Carlo simulations in %.2f seconds.\n',numTraj, time);

end

function [ time, loc ] = trajectory( L, edges, edgeJumps, startNode, startNodeInd, tmax, negateSteps )
% must have length of unit cell = 1

% check if start is valid
if nargin < 5 || isempty(startNodeInd)
    startNodeInd = 1;
end
% default: choose tmax so that we take 1000 steps
if nargin < 6 || isempty(tmax)
    avgRate = -mean(diag(L));
    tmax = 1000/avgRate;
end

if nargin < 7 || isempty(negateSteps)
    negateSteps = 0;
end

dim = size(edgeJumps,2);
curInd = startNodeInd;

% Set up vector of time points
avgRate = -mean(diag(L));
approxNumSteps = ceil(1.25*tmax*avgRate);
time = zeros(1,approxNumSteps);
node = zeros(1,approxNumSteps);

stepNum = 1;
t = 0;
while t < tmax
    % update time
    lambda = -L(curInd,curInd);
    dt = exprnd(1/lambda);
    t = t + dt;
    time(stepNum + 1) = t;

    % rates of leaving
    rates = L(curInd,:);
    rates(rates < 0) = 0;
    
    % update pos
    dxRand = lambda*rand(1);
    newInd = find(cumsum(rates) >= dxRand,1);
    
    node(stepNum) = newInd;
    
    curInd = newInd;
    stepNum = stepNum + 1;
end

node = node(1:stepNum - 1);
edges_traj = zeros(stepNum-1,2);
edges_traj(:,1) = [startNodeInd; node(1:end-1)'];
edges_traj(:,2) = node';

[~,edgeInds_traj] = ismember(edges_traj,edges,'rows');
time = time(1:stepNum);

edgeJumps_traj = edgeJumps(edgeInds_traj,:);

% Only use if m = 2
if negateSteps
    stepsToNegate = rand(stepNum-1,1) < .5;
    edgeJumps_traj(stepsToNegate,:) = -edgeJumps_traj(stepsToNegate,:);
end

loc = zeros(dim,stepNum);
loc(:,1) = startNode;
loc(:,2:end) = cumsum(edgeJumps_traj,1)' + startNode';

end


function [Deff, DeffVar, CI_high, CI_low] = getDeff(t,sd)

[N, T] = size(sd);

Ybar = mean(sd,2);
alpha = t./sum(t.^2);
slope = alpha*Ybar;

covYbar = zeros(T,T);

for i = 1:T
    for j = 1:T
        Yi = sd(:,i);
        Yj = sd(:,j);
        Yibar = mean(Yi);
        Yjbar = mean(Yj);
        
        covYbar(i,j) = (1/N)*(1/N*sum(Yi.*Yj)-Yibar*Yjbar);
    end
end

%varslope = alpha*covYbar*alpha';

Deff = slope/(2*d*t(end));
DeffVar = alpha*covYbar*alpha'/((2*d*t(end))^2);
CI_low = Deff - 1.96*sqrt(DeffVar);
CI_high = Deff + 1.96*sqrt(DeffVar);
end

%{
function [Deff, CI_high CI_low] = getDeff(t,sd)

[N, T] = size(sd);

Ybar = mean(sd,2);
alpha = t./sum(t.^2);
slope = alpha*Ybar;

for i = 1:T
    for j = 1:T
        Yi = sd(:,i);
        Yj = sd(:,j);
        Yibar = mean(Yi);
        Yjbar = mean(Yj);
        
        covYbar(i,j) = (1/N)*(1/N*sum(Yi.*Yj)-Yibar*Yjbar);
    end
end

varslope = alpha*covYbar*alpha';

Deff = slope/(2*d*t(end));
CI_high = Deff+1.96*sqrt(varslope)/(2*d*t(end))
end
%}