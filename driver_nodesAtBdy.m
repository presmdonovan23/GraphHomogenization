% runs square obstruction with

saveOn = 1;
plotOn = 0;
rho = .5;
mVals = 2.^[2:9];
geometry = 'squareBdySlow';
delta = 2;
dim = 2;
D0 = 1;

rate = [];

%% no diag jumps
diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
i = 1;
for m = mVals
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    rateCoeffs.alpha = 0;
    rateCoeffs.K1 = [];
    rateCoeffs.K2 = [];
    rateCoeffs.delta = delta;
    
    ctr = .75 - .5*h;
    %ctr = .5 - .5*h;
    rhoCorrected = rho - .5*h;

    ghParams(i) = GraphHomogParams_lattice(dim,geometry,D0,rhoCorrected,m,rate,rateCoeffs,diagJumps,ctr);
    ghInput(i) = GraphHomogInput(ghParams(i));
    ghParams(i).rho = rho;

    results_homog(i) = getDeff_homog(ghInput(i));
    results_mc = [];
    if plotOn
        plotGraph_lattice(ghParams(i),ghInput(i));
    end
    i = i+1;
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

%% diag jumps
diagJumps = 1; % diagJumps = 2 => corrected diagonal jumps

i = 1;
for m = mVals
    h = 1/m;
    
    rateCoeffs.dist = .9999*h;
    rateCoeffs.alpha = 0;
    rateCoeffs.K1 = [];
    rateCoeffs.K2 = [];
    rateCoeffs.delta = delta;
    
    ctr = .75 - .5*h;
    %ctr = .5 - .5*h;
    rhoCorrected = rho - .5*h;

    ghParams(i) = GraphHomogParams_lattice(dim,geometry,D0,rhoCorrected,m,rate,rateCoeffs,diagJumps,ctr);
    ghInput(i) = GraphHomogInput(ghParams(i));
    ghParams(i).rho = rho;

    results_homog(i) = getDeff_homog(ghInput(i));
    results_mc = [];
    if plotOn
        plotGraph_lattice(ghParams(i),ghInput(i));
    end
    i = i+1;
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
