function fh = plotGraph_lattice(ghParams,ghInput,plotEdges,plotObs)

if nargin < 3 || isempty(plotEdges)
    plotEdges = 1;
end
if nargin < 4 || isempty(plotObs)
    plotObs = 1;
end

nodes = ghInput.nodes;
edges = ghInput.edges;
edgeRates = ghInput.edgeRates;
edgeJumps = ghInput.edgeJumps;

dim = ghParams.dim;
geometry = ghParams.geometry;
obRad = ghParams.R;
maxHeadSize = max(sqrt(sum(edgeJumps.^2,2)));

fh = figure;
hold on

if dim == 2
    plot(nodes(:,1),nodes(:,2),'b.','markersize',35);
    
    if plotEdges
        quiver( nodes(edges(:,1),1),nodes(edges(:,1),2),...
                edgeJumps(:,1),edgeJumps(:,2),...
                'color','black','autoscale','off','MaxHeadSize',maxHeadSize);
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
elseif dim == 3
    plot3(nodes(:,1),nodes(:,2),nodes(:,3),'b.','markersize',35);
    
    if plotEdges
        quiver3(nodes(edges(:,1),1),nodes(edges(:,1),2),nodes(edges(:,1),3),...
                edgeJumps(:,1),edgeJumps(:,2),edgeJumps(:,3),...
                'color','black','autoscale','off','MaxHeadSize',maxHeadSize);
    end
    
    if plotObs
        if strcmpi(geometry,'circle')
            center = .5;
            
            [X,Y,Z] = sphere(100);
            X = X*obRad + center;
            Y = Y*obRad + center;
            Z = Z*obRad + center;
            surf(X,Y,Z,'facecolor','r','FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
            
        elseif strcmpi(geometry,'square')
            center = [.5 .5 .5];
            sideLength = obRad*2;
            myDrawCube(center,sideLength,fh);
        end
    end
end

% only label jumps in 2d
if dim == 2
    for i = 1:size(edges,1)
        edge = edges(i,:);
        rate = edgeRates(i);
        jump = edgeJumps(i,:);

        % plot rates as text
        textloc = nodes(edge(1),:) + .5*jump;
        if jump(1) > 10e-6 %edge goes left to right
            textloc(2) = textloc(2) + norm(jump)/10;

            h = text(textloc(1),textloc(2),num2str(rate,'%.2f'),...
                'horizontalalignment','center');

        elseif jump(1) < -10e-6 %edge goes right to left
            textloc(2) = textloc(2) - norm(jump)/10;

            h = text(textloc(1),textloc(2),num2str(rate,'%.2f'),...
                'horizontalalignment','center');
            set(h, 'rotation', 180)
        elseif jump(2) > 10e-6 %edge goes bottom to top
            textloc(1) = textloc(1) - norm(jump)/10;

            h = text(textloc(1),textloc(2),num2str(rate,'%.2f'),...
                'horizontalalignment','center');
            set(h, 'rotation', 90)
        elseif jump(2) < -10e-6 %edge goes top to bottom
            textloc(1) = textloc(1) + norm(jump)/10;

            h = text(textloc(1),textloc(2),num2str(rate,'%.2f'),...
                'horizontalalignment','center');
            set(h, 'rotation', 270)
        end
    end
end

end