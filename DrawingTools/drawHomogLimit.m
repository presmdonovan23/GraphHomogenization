function fh = drawHomogLimit(rho,geometryName,m,epsInvVals,drawObs,saveOn)
% ---Calling syntax---
%
% rho = .45;
% geometryName = 'circle';
% m = 16;
% epsInvVals = [1 2 4];%[1 2 4 8];
% drawObs = 0;
% 
% drawHomogLimit(rho,geometryName,m,epsInvVals,drawObs,saveOn)

if nargin < 4 || isempty(epsInvVals)
    epsInvVals = [1, 2, 4];
end
if nargin < 5 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 6 || isempty(saveOn)
    saveOn = 0;
end

h = 1/m;

ctr = [];
drawEdges = 1;
drawRates = 1;

fh = figure;
idx = 1;
for epsInv = epsInvVals
    subplot(1,length(epsInvVals),idx);

    for cellNumX = 1:epsInv
        for cellNumY = 1:epsInv
            loc = [cellNumX, cellNumY] - 1;
            
            fh = drawCell_lattice(rho,m,ctr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh);
        end
    end
    axis square
    axis([-h*epsInv epsInv+h*epsInv -h*epsInv epsInv+h*epsInv])
    axis off
   
    idx = idx + 1;
end

if saveOn
    filename = ['homogLimit_' geometryName '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end