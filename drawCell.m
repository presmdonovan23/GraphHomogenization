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
    
    ctr = ghParams.ctr;
    
    if strcmpi(geometry,'circle')
        obRad = rho/2;
        cornerX = loc(1) + ctr - obRad;
        cornerY = loc(2) + ctr - obRad;
        pos = [cornerX cornerY 2*obRad 2*obRad];
        rectangle(  'Position', pos,...
                    'facecolor','red',...
                    'Curvature',[1 1]);
        
    elseif strcmpi(geometry,'square')
        s = rho;

        cornerX = loc(1) + ctr - .5*s;
        cornerY = loc(2) + ctr - .5*s;

        rectangle('Position',[cornerX cornerY s s],'edgecolor','red','facecolor','red');
    end
end

if drawEdges
    
    if ghParams.m == 2
        warning('Hardcoded modification for m = 2.');
        nEdgeJumps = size(ghInput.edgeJumps,1);
        [~,inds] = unique([ghInput.nodes(ghInput.edges(:,1),:),ghInput.edgeJumps],'rows');
        negInds = setdiff((1:nEdgeJumps),inds);
        ghInput.edgeJumps(negInds,:) = -ghInput.edgeJumps(negInds,:);
        
    end
    for i = 1:size(ghInput.edgeJumps)
        edgeStart = ghInput.nodes(ghInput.edges(i,1),:);
        nu = ghInput.edgeJumps(i,:);

        myDrawLine(loc + edgeStart,loc + edgeStart + nu,4*3/ghParams.m,[],fh);
    end
end

nodesize = 2;
for i = 1:size(ghInput.nodes,1)
    
    cornerX = loc(1) + ghInput.nodes(i,1) - .5*nodesize*nodeRad;
    cornerY = loc(2) + ghInput.nodes(i,2) - .5*nodesize*nodeRad;
    pos = [cornerX cornerY nodeRad*nodesize nodeRad*nodesize];
    rectangle(  'Position', pos,...
                'facecolor','blue',...
                'edgecolor','blue',...
                'Curvature',[1 1]);

end

rectangle('Position',[loc(1)+0 loc(2)+0 1 1],'linewidth',1);

end