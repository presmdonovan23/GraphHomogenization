
saveOn = 0;

rhoVals = .5;
mVals = 2.^[2 3 4 5 6];
%deltaVals = [1 .9 .5 .25 .1 .05 .01 .001 .0001]; % squareSlowdown slowdown coeff

dim = 2;
geometry = 'square';

diagJumps = 0;
D0 = 1;
alpha = 0;

K1 = 25;
K2 = 10;

startNodeInd = 1;
numTraj = 0;
plotOn = 0;
rate = []; % field is current obsolete

idx = 1;
for m = mVals
    for rho = rhoVals
        %for delta = deltaVals

            rateCoeffs.alpha = alpha;
            rateCoeffs.K1 = 25;
            rateCoeffs.K2 = 10;
            %rateCoeffs.delta = delta;

            h = 1/m;

            ghParams(idx) = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
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
    if diagJumps
        dirName = [dirName '_diagJumps'];
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