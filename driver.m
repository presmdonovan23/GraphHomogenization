%% A simple driver
plotOn = 1;
saveOn = 0;
verbose = 1;

% geometry parameters
dim = 3;                    % dimension
m = 50;                      % mesh parameter
name = 'circle';            % obstruction type
obRad = .25;                % obstruction radius
obCtr = [.5 .5 .5];            % obstruction center
diagJumps = 0;              % diagonal jumps on/off
specialSetting = 'none';    % special setting (see README)
driftMult = 0;              % drift field strength
driftDecay = [];            % drift field decay
obSlowdownFctr = [];        % used for some special settings
bdyDist = [];               % used for some special settings

% monte carlo parameters
numTraj = 100;              % number of monte carlo simulations
startNodeInd = 1;           % inittial condition for monte carlo

% wrap up geometry parameters into an object
latticeGeo = LatticeGeometry(dim, m, name, obRad, ...
            obCtr, diagJumps, specialSetting, ...
            driftMult, driftDecay, obSlowdownFctr, bdyDist);

% convert geometry parameters into the necessary homogenization theory
% inputs
[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo, verbose);

% compute effective diffusivity
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,numTraj,startNodeInd,latticeGeo,verbose);
                
if saveOn
    saveResults(results);
end

if plotOn
    drawCell_lattice(latticeGeo);
end