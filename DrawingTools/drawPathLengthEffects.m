function fh = drawPathLengthEffects(rho,geometryName,mVals,drawObs,saveOn)
% ---Calling syntax---
%
% saveOn = 0;
% rho = .5;
% geometryName = 'square';
% mVals = [4 8 16 32];
% drawObs = 1; 
% drawPathLengthEffects(rho,geometryName,mVals,drawObs,saveOn)

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

ctr = [];
loc = [];
rate = [];
rateCoeffs = [];
drawEdges = 1;
drawRates = 0;

fh = figure;
idx = 1;
for m = mVals

    subplot(1,length(mVals),idx);
    
    fh = drawCell_lattice(rho,m,ctr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh);
    axis([-.5/mVals(1) 1+.5/mVals(1) -.5/mVals(1) 1+.5/mVals(1)])
   
    idx = idx + 1;
end

if saveOn
    filename = ['pathLengthEffects_' geometryName '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end