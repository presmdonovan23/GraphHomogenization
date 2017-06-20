
%myrate = @(x,j) rate_circle(x,j,h,D0,alpha,ctr,rho);
%myrate = @(x,j) rate_square_stick(x,h,D0,rho);
%myrate = @(x,j) rate_square_stick(x,h,D0,rho);
saveOn = 1;

rhoVals = [.5];
mVals = 4;%2.^[2:9];

dim = 2;
geometry = 'circle';
diagJumps = 0;
D0 = 1;
alpha = 0;
K1 = 25;
K2 = 10;

startNodeInd = 1;
numTraj = 0;
plotOn = 0;
rate = []; % field is current obsolete
rateCoeffs.alpha = alpha;
rateCoeffs.K1 = 25;
rateCoeffs.K2 = 10;

idx = 1;
for m = mVals
    for rho = rhoVals
        
        h = 1/m;
        
        ghParams(idx) = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
        ghInput = GraphHomogInput(ghParams(idx));
        
        results_homog(idx) = getDeff_homog(ghInput);
        results_mc(idx) = getDeff_MC( ghInput, numTraj, startNodeInd, plotOn );
        
        idx = idx+1;
    end
end

if saveOn
    dirname = 'Results';
    filename = sprintf('%s/results_%dd_%s',dirname,dim,geometry);
    if diagJumps
        filename = [filename '_diagJumps'];
    end
    filenameFull = myClockFilename(filename);
    save(filenameFull,'results_homog','results_mc','ghParams');
end