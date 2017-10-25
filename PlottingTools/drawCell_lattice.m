function [fh, nodes, edges, edgeJumps] = drawCell_lattice(rho,m,ctr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh)
% ---Calling syntax---
%
% saveOn = 0;
% rho = .45;
% m = 8;
% geometryName = 'circle';
% drawObs = 1;
% ctr = .5;
% drawSetting(rho,m,ctr,geometryName,drawObs,saveOn)

if nargin < 5 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 6 || isempty(drawEdges)
    drawEdges = 1;
end
if nargin < 7 || isempty(drawRates)
    drawEdges = 1;
end
if nargin < 8 || isempty(saveOn)
    saveOn = 0;
end
if nargin < 9 || isempty(loc)
    loc = [0, 0];
end
if nargin < 10 || isempty(fh)
    fh = figure;
end

dim = 2;

diagJumps = 0;
D0 = 1;
rateCoeffs = [];
specialSetting_m2 = 'none';

rate = [];

[~,nodes,edges,edgeRates,edgeJumps] = homogInputs_lattice(...
    dim,geometryName,D0,rho,m,rate,rateCoeffs,diagJumps,ctr,specialSetting_m2);
        
drawCell(loc,nodes,edges,edgeRates,edgeJumps,geometryName,m,rho,ctr,drawObs,drawEdges,drawRates,fh);

if saveOn
    filename = [geometryName '_rho' num2str(round(100*rho)) '_m' num2str(m)];
    mySaveFig(filename, fh, 'fig' );
    mySaveFig(filename, fh, 'png' );
    mySaveFig(filename, fh, 'eps' );
end

end