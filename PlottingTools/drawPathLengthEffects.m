function fh = drawPathLengthEffects(rho,geometryName,mVals,drawObs,obCtr,saveOn)
% ---Calling syntax---
%
% rho = .5;
% geometryName = 'square';
% mVals = [2 4 8 16];
% drawObs = 1; 
% obCtr = .75;
% saveOn = 0;
% drawPathLengthEffects(rho,geometryName,mVals,drawObs,obCtr,saveOn)

if nargin < 3 || isempty(mVals)
    mVals = [4 8 16];
end
if nargin < 4 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 5 || isempty(obCtr)
    obCtr = .5;
end
if nargin < 6 || isempty(saveOn)
    saveOn = 0;
end

loc = [];
drawEdges = 1;
drawRates = 0;

fh = figure;
idx = 1;
for m = mVals

    subplot(1,length(mVals),idx);
    
    fh = drawCell_lattice(rho,m,obCtr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh);
    axis([-.5/mVals(1) 1+.5/mVals(1) -.5/mVals(1) 1+.5/mVals(1)])
   
    idx = idx + 1;
end

if saveOn
    filename = ['pathLengthEffects_' geometryName '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end