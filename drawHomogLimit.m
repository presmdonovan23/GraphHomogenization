function fh = drawHomogLimit(rho,geometry,m,epsInvVals,drawObs,saveOn)

if nargin < 4 || isempty(epsInvVals)
    epsInvVals = [1, 2, 4];
end
if nargin < 5 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 6 || isempty(saveOn)
    saveOn = 0;
end

dim = 2;
diagJumps = 0;
D0 = 1;
h = 1/m;

alpha = 0;
rateCoeffs.alpha = alpha;

fh = figure;
idx = 1;
for epsInv = epsInvVals
    subplot(1,length(epsInvVals),idx);
    if strcmpi(geometry,'circle')
        rate = @(x,j) rate_circle(x,j,h,D0,alpha,rho);
    elseif strcmpi(geometry,'square')
        rate = @(x,j) rate_square(x,j,h,D0,alpha,rho);
    end
    ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
    ghInput = GraphHomogInput(ghParams);

    for cellNumX = 1:epsInv
        for cellNumY = 1:epsInv
            loc = [cellNumX, cellNumY] - 1;
            drawCell(loc,ghInput,ghParams,drawObs,fh);
        end
    end
    axis square
    axis([-h*epsInv epsInv+h*epsInv -h*epsInv epsInv+h*epsInv])
    axis off
   
    idx = idx + 1;
end

if saveOn
    filename = ['homogLimit_' geometry '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end