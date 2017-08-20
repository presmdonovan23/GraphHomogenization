
saveOn = 0;

rhoVals = .5;
mVals = 2.^[12];
%deltaVals = [1 .9 .5 .25 .1 .05 .01 .001 .0001]; % squareSlowdown slowdown coeff
ctr = .5;
dim = 2;
geometry = 'circle';

diagJumps = 0; % diagJumps = 2 => corrected diagonal jumps
D0 = 1;

alpha = 1;
K1 = 25;
K2 = 10;
delta = 1.5;

startNodeInd = 1;
numTraj = 0;
plotOn = 0;
rate = []; % field is current obsolete

clear ghParams ghInput results_homog results_mc
idx = 1;
for m = mVals
    for rho = rhoVals
        %for delta = deltaVals
            h = 1/m;
            
            rateCoeffs.alpha = alpha;
            rateCoeffs.K1 = 25;
            rateCoeffs.K2 = 10;
            rateCoeffs.delta = delta;

            ghParams(idx) = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps,ctr);
            ghInput(idx) = GraphHomogInput(ghParams(idx));

            results_homog(idx) = getDeff_homog(ghInput(idx));
            results_mc(idx) = getDeff_MC( ghInput(idx), numTraj, startNodeInd, plotOn );

            idx = idx+1;
        %end
    end
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