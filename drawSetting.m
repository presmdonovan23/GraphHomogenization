function fh = drawSetting(rho,m,geometry,drawObs,saveOn)

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

h = 1/m;

epsInv = 1;

if strcmpi(geometry,'circle')
    rate = @(x,j) rate_circle(x,j,h,D0,alpha,rho);
elseif strcmpi(geometry,'square')
    rate = @(x,j) rate_square(x,j,h,D0,alpha,rho);
end

ghParams = GraphHomogParams_lattice(dim,geometry,D0,rho,m,rate,rateCoeffs,diagJumps);
ghInput = GraphHomogInput(ghParams);

fh = figure;
hold on

%% draw initial things to get legend entries correct

%plot(.5,.5,'k.','markersize',20)
%plot(.5,.5,'k-','linewidth',1)
%plot(.5,.5,'r.','markersize',20)
%lh = legend('Node','Edge','Obstruction (no nodes)','location','eastoutside');
%set(lh,'fontsize',22);
%%
for cellNumX = 1:epsInv
    for cellNumY = 1:epsInv
        loc = [cellNumX, cellNumY] - 1;
        
        drawCell(loc,ghInput,ghParams,drawObs,fh);
    end
end
axis square
axis([-h*epsInv epsInv+h*epsInv -h*epsInv epsInv+h*epsInv])
axis off

if saveOn
    filename = ['setting_' geometry '_rho' num2str(round(100*rho))];
    mysavefig(filename, fh, 'fig' );
    mysavefig(filename, fh, 'png' );
end

end