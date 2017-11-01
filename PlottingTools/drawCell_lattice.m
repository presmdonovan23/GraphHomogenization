function [fh, nodes, edges, edgeJumps] = drawCell_lattice(latticeGeo,drawObs,drawEdges,drawRates,saveOn,loc,fh)
% ---Calling syntax---
%
% saveOn = 0;
% rho = .45;
% m = 8;
% geometryName = 'circle';
% drawObs = 1;
% ctr = .5;
% drawSetting(rho,m,ctr,geometryName,drawObs,saveOn)

if nargin < 2 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 3 || isempty(drawEdges)
    drawEdges = 1;
end
if nargin < 4 || isempty(drawRates)
    drawRates = 1;
end
if nargin < 5 || isempty(saveOn)
    saveOn = 0;
end
if nargin < 6 || isempty(loc)
    loc = zeros(1,latticeGeo.dim);
end
if nargin < 7 || isempty(fh)
    fh = figure;
end

[~,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(latticeGeo);
        
drawCell(latticeGeo,nodes,edges,edgeRates,edgeJumps,loc,drawObs,drawEdges,drawRates,fh);

if saveOn
    filename = [geometryName '_obRad' num2str(round(100*obRad)) '_m' num2str(m)];
    mySaveFig(filename, fh, 'fig' );
    mySaveFig(filename, fh, 'png' );
    mySaveFig(filename, fh, 'eps' );
end

end