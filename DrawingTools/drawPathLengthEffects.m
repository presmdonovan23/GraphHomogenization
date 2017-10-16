function fh = drawPathLengthEffects(rho,geometry,mVals,drawObs,saveOn)

if nargin < 3 || isempty(mVals)
    mVals = [4 8 16];
end
if nargin < 4 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 5 || isempty(saveOn)
    saveOn = 0;
end

dim = 2;

diagJumps = 0;
D0 = 1;
alpha = 0;
rateCoeffs.alpha = alpha;

fh = figure;
idx = 1;
for m = mVals

    h = 1/m;

    if strcmpi(geometry,'circle')
        rate = @(x,j) rate_circle(x,j,h,D0,alpha,rho);
    elseif strcmpi(geometry,'square')
        rate = @(x,j) rate_square(x,j,h,D0,alpha,rho);
    end
    ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
    ghInput = GraphHomogInput(ghParams);

    subplot(1,length(mVals),idx);
    
    drawCell([0 0],ghInput,ghParams,drawObs,fh);
        
    axis square
    axis([-.5/mVals(1) 1+.5/mVals(1) -.5/mVals(1) 1+.5/mVals(1)])
    axis off
   
    idx = idx + 1;
end

if saveOn
    filename = ['pathLengthEffects_' geometry '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end