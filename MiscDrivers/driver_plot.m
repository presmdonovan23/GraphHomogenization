%% A simple driver for plotting
plotOn = 1;
saveOn = 0;
verbose = 1;

% geometry parameters
dim = 2;                    % dimension
m = 5;                      % mesh parameter
name = 'square';            % obstruction type
obRad = .25;                % obstruction radius
obCtr = [.75 .75];            % obstruction center
diagJumps = 1;              % diagonal jumps on/off
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

if plotOn
    %drawCell_lattice(latticeGeo);
    
    drawEdges = 1;
    drawRates = 0;
    loc = zeros(1,dim);
    fh = figure;
    drawCell_lattice(latticeGeo,drawObs,drawEdges,drawRates,saveOn,loc,fh)
end