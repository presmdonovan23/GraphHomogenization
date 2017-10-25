function fh = drawHomogLimit(rho,geometryName,m,epsVals,drawObs,saveOn)
% Currently there is a minor issue where the node sizes will not shrink
% with epsilon. This can be easily implemented by allowing the node size to
% be passed to drawCell_lattice.
% ---Calling syntax---
%
% rho = .45;
% geometryName = 'circle';
% m = 16;
% epsVals = [1 1/2 1/4];
% drawObs = 0;
% saveOn = 0;
% drawHomogLimit(rho,geometryName,m,epsVals,drawObs,saveOn);

if nargin < 4 || isempty(epsVals)
    epsVals = [1, 1/2, 1/4];
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
drawRates = 0;

fh = figure;
idx = 1;
for i = 1:length(epsVals)
    epsInv = 1/epsVals(i);
    subplot(1,length(epsVals),idx);

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

% nodes look better if their size is scaled by epsilon
for i = 1:length(fh.Children)
    
    eps = epsVals(1 + length(fh.Children) - i);
    
    markerChildren = arrayfun(@(c) isa(c,'matlab.graphics.chart.primitive.Line'),fh.Children(i).Children);
    mSize = 1.25*fh.Children(i).Children(find(markerChildren,1)).MarkerSize*eps; % get marker size of first valid child
    set(fh.Children(i).Children(markerChildren),'MarkerSize',mSize);
    
end
if saveOn
    filename = ['homogLimit_' geometryName '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end