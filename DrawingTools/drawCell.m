function fh = drawCell(loc,nodes,edges,edgeRates,edgeJumps,geometryName,m,rho,ctr,drawObs,drawEdges,drawRates,fh)

dim = size(nodes,2);
if nargin < 1 || isempty(loc)
    loc = zeros(1,dim);
end
if nargin < 9 || isempty(ctr)
    ctr = .5;
end
if nargin < 10 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 11 || isempty(drawEdges)
    drawEdges = 1;
end
if nargin < 12 || isempty(drawRates)
    drawRates = 0;
end
if nargin < 13 || isempty(fh)
    fh = figure;
end

figure(fh);
hold on;

h = 1/m;
nodeRad = 175*h;
obRad = rho/2;
if drawObs
    
    if dim == 2
        if contains(geometryName,'circle')
            
            cornerX = loc(1) + ctr - obRad;
            cornerY = loc(2) + ctr - obRad;
            pos = [cornerX cornerY 2*obRad 2*obRad];
            rectangle(  'Position', pos,...
                        'facecolor','red',...
                        'Curvature',[1 1]);

        elseif contains(geometryName,'square')
            s = rho;

            cornerX = loc(1) + ctr - .5*s;
            cornerY = loc(2) + ctr - .5*s;

            rectangle('Position',[cornerX cornerY s s],'edgecolor','red','facecolor','red');
        end
    elseif dim == 3
        if contains(geometryName,'circle')
            center = loc + .5;
            
            [X,Y,Z] = sphere(100);
            X = X*obRad + center(1);
            Y = Y*obRad + center(2);
            Z = Z*obRad + center(3);
            surf(X,Y,Z,'facecolor','r','FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
            
        elseif contains(geometryName,'square')
            center = loc + [.5 .5 .5];
            sideLength = obRad*2;
            myDrawCube(center,sideLength,fh);
        end
    end
end

if drawEdges
    
    if m == 2
        warning('reminder: hard-coded modification for m = 2.');
        nEdgeJumps = size(edgeJumps,1);
        [~,inds] = unique([nodes(edges(:,1),:),edgeJumps],'rows');
        negInds = setdiff((1:nEdgeJumps),inds);
        edgeJumps(negInds,:) = -edgeJumps(negInds,:);
        
    end
    
    if drawEdges == 1 % draw edges w/o arrows
        for i = 1:size(edgeJumps)
            edgeStart = nodes(edges(i,1),:);
            nu = edgeJumps(i,:);
        
            myDrawLine(loc + edgeStart,loc + edgeStart + nu,4*3/m,[],fh);
        end
    elseif drawEdges == 2 % draw edges w/ arrows
        maxHeadSize = .5*max(sqrt(sum(edgeJumps.^2,2)));
        if dim == 2
            quiver( loc(1) + nodes(edges(:,1),1),loc(2) + nodes(edges(:,1),2),...
                    edgeJumps(:,1),edgeJumps(:,2),...
                    'color','black','autoscale','off','MaxHeadSize',maxHeadSize);
        elseif dim == 3
            quiver3( loc(1) + nodes(edges(:,1),1),loc(2) + nodes(edges(:,1),2),loc(3) + nodes(edges(:,1),3),...
                    edgeJumps(:,1),edgeJumps(:,2),edgeJumps(:,3),...
                    'color','black','autoscale','off','MaxHeadSize',maxHeadSize);
        end
    end
    
    
end

%% plot rates
if drawRates && dim == 2
    for i = 1:size(edges,1)
        edge = edges(i,:);
        rate = edgeRates(i);
        jump = edgeJumps(i,:);

        % plot rates as text
        textloc = loc + nodes(edge(1),:) + .5*jump;
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

%% plot nodes
if dim == 2
    plot(nodes(:,1),nodes(:,2),'b.','markersize',nodeRad);
elseif dim == 3
    plot3(loc(1) + nodes(:,1),loc(2) + nodes(:,2),loc(3) + nodes(:,3),'b.','markersize',35);
end

%% plot boundary
if dim == 2
    rectangle('Position',[loc(1)+0 loc(2)+0 1 1],'linewidth',1);
end

axis square
axis off


end