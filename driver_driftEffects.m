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

startNodeInd = 1;
numTraj = 0;
plotOn = 1;
rate = []; % field is current obsolete

rateCoeffs.dist = .9999*h;
rateCoeffs.alpha = 0;
rateCoeffs.K1 = [];
rateCoeffs.K2 = [];
%% regular
geometry = 'square';

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

results_homog = getDeff_homog(ghInput);
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );
results_regular = results_homog;
if plotOn
    plotGraph_lattice(ghParams,ghInput);
end
if saveOn
    fileName = sprintf('results_%s',myClock(6));
    
    dirName = sprintf('Results_%dd_%s',dim,geometry);
    if diagJumps == 1
        dirName = [dirName '_diagJumps'];
    elseif diagJumps == 2
        dirName = [dirName '_diagJumpsCorrected'];
    end
    
    if ~exist(dirName,'dir')
        mkdir(dirName);
    end
    %filename = sprintf('%s/results_%dd_%s',dirname,dim,geometry);
    %filenameFull = myClockFilename(filename);
    filenameFull = [dirName '/' fileName];
    fprintf('Filename: %s\n',filenameFull);
    save(filenameFull,'results_homog','results_mc','ghParams','ghInput');
end

%% bonding
geometry = 'squareBonding';
delta = .5;
rateCoeffs.delta = delta;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

results_homog = getDeff_homog(ghInput);
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );
results_bonding = results_homog;
if plotOn
    plotGraph_lattice(ghParams,ghInput);
end
if saveOn
    fileName = sprintf('results_%s',myClock(6));
    
    dirName = sprintf('Results_%dd_%s',dim,geometry);
    if diagJumps == 1
        dirName = [dirName '_diagJumps'];
    elseif diagJumps == 2
        dirName = [dirName '_diagJumpsCorrected'];
    end
    
    if ~exist(dirName,'dir')
        mkdir(dirName);
    end
    %filename = sprintf('%s/results_%dd_%s',dirname,dim,geometry);
    %filenameFull = myClockFilename(filename);
    filenameFull = [dirName '/' fileName];
    fprintf('Filename: %s\n',filenameFull);
    save(filenameFull,'results_homog','results_mc','ghParams','ghInput');
end

%% boundary attraction
geometry = 'squareBdyAttract';
delta = 2;
rateCoeffs.delta = delta;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

results_homog = getDeff_homog(ghInput);
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );
results_bdyAttract = results_homog;
if plotOn
    plotGraph_lattice(ghParams,ghInput);
end
if saveOn
    fileName = sprintf('results_%s',myClock(6));
    
    dirName = sprintf('Results_%dd_%s',dim,geometry);
    if diagJumps == 1
        dirName = [dirName '_diagJumps'];
    elseif diagJumps == 2
        dirName = [dirName '_diagJumpsCorrected'];
    end
    
    if ~exist(dirName,'dir')
        mkdir(dirName);
    end
    %filename = sprintf('%s/results_%dd_%s',dirname,dim,geometry);
    %filenameFull = myClockFilename(filename);
    filenameFull = [dirName '/' fileName];
    fprintf('Filename: %s\n',filenameFull);
    save(filenameFull,'results_homog','results_mc','ghParams','ghInput');
end

%% boundary repel
geometry = 'squareBdyRepel';
delta = 2;
rateCoeffs.delta = delta;

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

results_homog = getDeff_homog(ghInput);
results_mc = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );
results_bdyRepel = results_homog;
if plotOn
    plotGraph_lattice(ghParams,ghInput);
end
if saveOn
    fileName = sprintf('results_%s',myClock(6));
    
    dirName = sprintf('Results_%dd_%s',dim,geometry);
    if diagJumps == 1
        dirName = [dirName '_diagJumps'];
    elseif diagJumps == 2
        dirName = [dirName '_diagJumpsCorrected'];
    end
    
    if ~exist(dirName,'dir')
        mkdir(dirName);
    end
    %filename = sprintf('%s/results_%dd_%s',dirname,dim,geometry);
    %filenameFull = myClockFilename(filename);
    filenameFull = [dirName '/' fileName];
    fprintf('Filename: %s\n',filenameFull);
    save(filenameFull,'results_homog','results_mc','ghParams','ghInput');
end

%% plot results
if plotOn
    fh = figure;
    bh = bar([results_regular.Deff,results_bonding.Deff,results_bdyRepel.Deff,results_bdyAttract.Deff]);
    text(bh.XData,bh.YData',num2str(bh.YData','%0.2f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ca = gca;
    ca.XTickLabel = {'Simple','Bonding','Repulsion','Attraction'};
    ylabel('D_e')
    ca = gca;
    ca.FontSize = 18;
    
    %dirName = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';
    %filename = 'driftEffectsDeff';
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'fig')
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'png')
    %mySaveFig([dirName 'driftEffectsDeff'],fh,'eps')

end

