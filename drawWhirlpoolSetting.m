function fh = drawWhirlpoolSetting( ghParams, ghInput, plotEdges, plotObs, plotLambda )

if nargin < 3 || isempty(plotEdges)
    plotEdges = 1;
end
if nargin < 4 || isempty(plotObs)
    plotObs = 0;
end
if nargin < 5 || isempty(plotLambda)
    plotLambda = 0;
end

nodes = ghInput.nodes;
edges = ghInput.edges;
edgeRates = ghInput.edgeRates;
edgeJumps = ghInput.edgeJumps;

geometry = ghParams.geometry;
obRad = ghParams.R;
maxHeadSize = .35*max(sqrt(sum(edgeJumps.^2,2)));

fh = figure;
hold on

plot(nodes(:,1),nodes(:,2),'b.','markersize',50);

for i = 1:size(nodes,1)
    s = sprintf('y_%d',i);
    text(nodes(i,1) - .09,nodes(i,2) - .06,s,'fontsize',20);
end
if plotEdges
    quiver( nodes(edges(:,1),1),nodes(edges(:,1),2),...
            edgeJumps(:,1),edgeJumps(:,2),...
            'color','black',...
            'autoscale','off',...
            'MaxHeadSize',maxHeadSize,...
            'linewidth',3);
end
if plotObs

    cornerX = .5 - obRad;
    cornerY = .5 - obRad;
    pos = [cornerX cornerY 2*obRad 2*obRad];

    if strcmpi(geometry,'circle')
        curvature = [1 1];
    elseif strcmpi(geometry,'square')
        curvature = [0 0];
    end

    rectangle(  'Position', pos,...
                'facecolor','red',...
                'Curvature',curvature);

end


for i = 1:size(edges,1)
    edge = edges(i,:);
    rate = edgeRates(i);
    jump = edgeJumps(i,:);

    % plot rates as text
    textloc = nodes(edge(1),:) + .5*jump;
    %rateStr = sprintf('%d',round(rate));%num2str(rate,'%.2f')
    if plotLambda
        rateStr = '$\bar{\lambda}$';
    else
        rateStr = '';
    end
    if jump(1) > 10e-6 %edge goes left to right
        textloc(2) = textloc(2) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');

    elseif jump(1) < -10e-6 %edge goes right to left
        textloc(2) = textloc(2) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');

        %set(h, 'rotation', 180)
    elseif jump(2) > 10e-6 %edge goes bottom to top
        textloc(1) = textloc(1) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');

        %set(h, 'rotation', 90)
    elseif jump(2) < -10e-6 %edge goes top to bottom
        textloc(1) = textloc(1) + norm(jump)/10;

        h = text(textloc(1),textloc(2),rateStr,...
            'horizontalalignment','center',...
            'fontsize',30,...
            'interpreter','latex');

        %set(h, 'rotation', 270)
    end
end

lh = legend('Nodes','Edges','location','southeast');
set(lh,'fontsize',18);
axis square
axis off
rectangle('position',[0 0 1 1],'linewidth',1)
%mySaveFig('detailedBalance',fh);

end

