function fh = drawCell(loc,ghInput,ghParams,drawObs,fh, drawEdges)

if nargin < 4 || isempty(drawObs)
    drawObs = 0;
end
if nargin < 5 || isempty(fh)
    fh = figure;
end
if nargin < 6 || isempty(drawEdges)
    drawEdges = 1;
end

nodeRad = 1.75*ghParams.h/10;

geometry = ghParams.geometry;
rho = ghParams.rho;

if drawObs
    if strcmpi(geometry,'circle')
        obRad = rho/2;
        cornerX = loc(1) + .5 - obRad;
        cornerY = loc(2) + .5 - obRad;
        pos = [cornerX cornerY 2*obRad 2*obRad];
        rectangle(  'Position', pos,...
                    'facecolor','red',...
                    'Curvature',[1 1]);
        
    elseif strcmpi(geometry,'square')
        s = rho;

        cornerX = loc(1) + .5*(1-s);
        cornerY = loc(2) + .5*(1-s);

        rectangle('Position',[cornerX cornerY s s],'edgecolor','red','facecolor','red');
    end
end

for i = 1:size(ghInput.nodes,1)
    
    cornerX = loc(1) + ghInput.nodes(i,1) - nodeRad;
    cornerY = loc(2) + ghInput.nodes(i,2) - nodeRad;
    pos = [cornerX cornerY 2*nodeRad 2*nodeRad];
    rectangle(  'Position', pos,...
                'facecolor','black',...
                'Curvature',[1 1]);

end

if drawEdges
    for i = 1:size(ghInput.edgeJumps)
        edgeStart = ghInput.nodes(ghInput.edges(i,1),:);
        nu = ghInput.edgeJumps(i,:);

        drawLine(loc + edgeStart,loc + edgeStart + nu,[],[],fh);
    end
end

rectangle('Position',[loc(1)+0 loc(2)+0 1 1],'linewidth',3);

end