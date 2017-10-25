% runs square obstruction with:
%   1) regular rates
%   2) bonding at boundary
%   3) attraction at boundary
%   4) repulsion at boundary

saveOn = 0;
rho = .5;
m = 8;
h = 1/m;
dim = 2;

diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
D0 = 1;
ctr = .5;
specialSetting_m2 = 'none';

startNodeInd = 1;
numTraj = 0;
plotOn = 1;
rate = []; % field is current obsolete
drawObs = 1;
drawEdges = 1;
drawRates = 1;

rateCoeffs.dist = .9999*h;
rateCoeffs.alpha = 0;
rateCoeffs.K1 = [];
rateCoeffs.K2 = [];
geometry.dim = dim;
geometry.m = m;
geometry.rho = rho;
geometry.diagJumps = diagJumps;
geometry.ctr = ctr;
geometry.specialSetting_m2 = specialSetting_m2;

%% regular
geometry.name = 'square';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,rho,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{1} = saveResults(results);
end

%% bonding
rateCoeffs.delta = .5;
geometry.name = 'squareBonding';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn

    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,rho,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{2} = saveResults(results);
end

%% boundary attraction
rateCoeffs.delta = 2;
geometry.name = 'squareBdyAttract';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,rho,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{3} = saveResults(results);
end

%% boundary repel
rateCoeffs.delta = 2;
geometry.name = 'squareBdyRepel';
geometry.rateCoeffs = rateCoeffs;

[L,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(dim,geometry.name,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
results = effDiff(L,nodes,edges,edgeRates,edgeJumps,[],[],geometry);

if plotOn
    
    drawEdges = 1;
    drawObs = 1;
    drawRates = 0;
    
    drawCell([],nodes,edges,edgeRates,edgeJumps,geometry.name,m,rho,ctr,drawObs,drawEdges,drawRates,[]);
    
end
if saveOn
    filename{4} = saveResults(results);
end

%% plot results
if plotOn
    load(filename{1});
    results_regular = results;
    load(filename{2});
    results_bonding = results;
    load(filename{3});
    results_bdyAttract = results;
    load(filename{4});
    results_bdyRepel = results;
    
    fh = figure;
    bh = bar([results_regular.Deff,results_bonding.Deff,results_bdyRepel.Deff,results_bdyAttract.Deff]);
    text(bh.XData,bh.YData',num2str(bh.YData','%0.2f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ca = gca;
    ca.XTickLabel = {'Neutral','Bonding','Repulsion','Attraction'};
    ylabel('D_e')
    ca = gca;
    ca.FontSize = 18;
    
    %dirName = '';
    %filename = 'driftEffectsDeff';
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'fig')
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'png')
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'eps')

end

