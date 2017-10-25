function fh = drawDriftField(rho,m,geometryName,rateCoeffs,drawObs,saveOn)
% ---Calling syntax---
%
% rho = .45;
% m = 8;
% geometryName = 'circle';
% rateCoeffs.alpha = 1;
% rateCoeffs.K1 = 25;
% rateCoeffs.K2 = 10;
% drawObs = 1;
% saveOn = 0;
% drawDriftField(rho,m,geometryName,rateCoeffs,drawObs,saveOn);

if nargin < 5 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 6 || isempty(saveOn)
    saveOn = 0;
end

h = 1/m;
ctr = .5;

fh = figure;
hold on

loc = [0,0];
drawEdges = 0;
drawRates = 0;

[fh, nodes] = drawCell_lattice(rho,m,ctr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh);
              
R = rho/2;
K1 = rateCoeffs.K1;
K2 = rateCoeffs.K2;

d = drift(nodes,ctr,R,K1,K2);
quiver(nodes(:,1),nodes(:,2),...
        nodes(:,1) + d(:,1),...
        nodes(:,2) + d(:,2),...
        'linewidth',3);

axis square
axis([-h 1+h -h 1+h])
axis off

if saveOn
    filename = ['driftField_' geometryName '_rho' num2str(round(100*rho)) '_alpha' num2str(alpha)];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end

function val = drift(x,ctr,R,K1,K2)
% only makes sense for circular obstuctions

relX = x - ctr;
relNorm = sqrt(sum(relX.^2,2));

dist2ob = relNorm - R;
xNormalized = relX./relNorm;

val = K1*xNormalized.*exp(-K2*dist2ob);
val(~isfinite(val)) = 0;

end