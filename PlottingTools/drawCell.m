function fh = drawCell(latticeGeo,nodes,edges,edgeRates,edgeJumps,loc,drawObs,drawEdges,drawRates,fh)

dim = latticeGeo.dim;

if nargin < 6 || isempty(loc)
    loc = zeros(1,dim);
end
if nargin < 7 || isempty(drawObs)
    drawObs = 1;
end
if nargin < 8 || isempty(drawEdges)
    drawEdges = 1;
end
if nargin < 9 || isempty(drawRates)
    drawRates = 0;
end
if nargin < 10 || isempty(fh)
    fh = figure;
end

m = latticeGeo.m;
obCtr = latticeGeo.obCtr;
obRad = latticeGeo.obRad;
name = latticeGeo.name;
nodeRad = 175*latticeGeo.h;

figure(fh);
hold on;

if drawObs
    
    if dim == 2
        if strcmpi(name,'circle')
            
            cornerX = loc(1) + obCtr(1) - obRad;
            cornerY = loc(2) + obCtr(2) - obRad;
            pos = [cornerX cornerY 2*obRad 2*obRad];
            rectangle(  'Position', pos,...
                        'facecolor','red',...
                        'Curvature',[1 1]);

        elseif strcmpi(name,'square')
            sideLen = latticeGeo.sideLen;

            cornerX = loc(1) + obCtr(1) - .5*sideLen;
            cornerY = loc(2) + obCtr(2) - .5*sideLen;

            rectangle('Position',[cornerX cornerY sideLen sideLen],'edgecolor','red','facecolor','red');
        end
    elseif dim == 3
        if strcmpi(name,'circle')
            center = loc + obCtr;
            
            [X,Y,Z] = sphere(100);
            X = X*obRad + center(1);
            Y = Y*obRad + center(2);
            Z = Z*obRad + center(3);
            surf(X,Y,Z,'facecolor','r','FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
            
        elseif strcmpi(name,'square')
            center = loc + obCtr;
            sideLength = obRad*2;
            myDrawCube(center,sideLength,fh);
        end
    end
end

if drawEdges
    
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
    plot(loc(1) + nodes(:,1),loc(2) + nodes(:,2),'b.','markersize',nodeRad);
elseif dim == 3
    plot3(loc(1) + nodes(:,1),loc(2) + nodes(:,2),loc(3) + nodes(:,3),'b.','markersize',nodeRad);
end

%% plot boundary
if dim == 2
    rectangle('Position',[loc(1)+0 loc(2)+0 1 1],'linewidth',1);
end

axis square
axis off


end

function fh = myDrawCube( center, sideLength, fh )
% code modified from a post on mathworks forums
if nargin < 3 || isempty(fh)
    fh = figure;
else
    figure(fh);
end

x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*sideLength+center(1);
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*sideLength+center(2);
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*sideLength+center(3);

for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),'r','FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none');
    set(h,'edgecolor','k')
end

end

function fh = myDrawLine( start,finish,width,color,fh )
%DRAWLINE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5 || isempty(fh)
    fh = figure;
else
    figure(fh);
end

if nargin < 3 || isempty(width)
    width = 1;
end

if nargin < 4 || isempty(color)
    color = 'black';
end

line([start(1) finish(1)],[start(2) finish(2)],...
            'linewidth',width,'color',color);

end