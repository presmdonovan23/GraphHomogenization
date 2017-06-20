function fh = drawDriftField(rho,m,geometry,alpha,drawObs,saveOn)

if nargin < 5 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 6 || isempty(saveOn)
    saveOn = 0;
end


dim = 2;

diagJumps = 0;
D0 = 1;

rateCoeffs.alpha = alpha;

h = 1/m;

if strcmpi(geometry,'circle')
    rate = @(x,j) rate_circle(x,j,h,D0,alpha,rho);
elseif strcmpi(geometry,'square')
    rate = @(x,j) rate_square(x,j,h,D0,alpha,rho);
end

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

fh = figure;
hold on

%%

loc = [0,0];
drawEdges = 0;
drawCell(loc,ghInput,ghParams,drawObs,fh,drawEdges);

nodes = ghInput.nodes;
drifts = drift_circle(nodes,[1,2],rho);
quiver(nodes(:,1),nodes(:,2),...
        nodes(:,1) + drifts(:,1),...
        nodes(:,2) + drifts(:,2),...
        'linewidth',3);
%for i = 1:size(nodes,1)
%    node = ghInput.nodes(i,:);
%    f_d = drift_circle(node,[1,2],rho);
%    quiver(node(1),node(2),node(1) + f_d(1),node(2) + f_d(2));
%end

axis square
axis([-h 1+h -h 1+h])
axis off

if saveOn
    filename = ['driftField_' geometry '_rho' num2str(round(100*rho)) '_alpha' num2str(alpha)];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end