function [nodes, isFree] = getNodes_lattice( R, m, dim, geometry, ctr )
% assumes nodes live on lattice
% assumes cell has length 1
tic

if nargin < 4
    geometry = 'circle';
end
if nargin < 5
    ctr = .5;
end

sz = [m*ones(1,dim) dim];
nodes = zeros(sz);
% store all possible sites on Z+
for i = 1:dim
    sz = ones(1,dim);
    sz(i) = m;
    v = zeros(sz);
    v(:) = 1:m;
    sz = m*ones(1,dim);
    sz(i) = 1;
    if dim == 2
        nodes(:,:,i) = repmat(v,sz);
    elseif dim == 3
        nodes(:,:,:,i) = repmat(v,sz);
    end
end
% shift/scale sites to live on [0,1]^d
nodes = (nodes-.5)*(1/m);

if R == 0
    mDim = m*ones(1,dim);
    isFree = true(mDim);
else
    if strcmpi(geometry,'circle')
        dist2ctr2 = sum((nodes-ctr).^2,dim+1);
        dist2ctr2 = round(dist2ctr2,m); % eliminate numerical error that could make obstructed sites non-symmetrical
        isFree = dist2ctr2 > R^2;
    elseif strcmpi(geometry,'square') % ** could be some numerical error for small m where sites are on boundary of obstructed region
        if dim == 2
            isFree = ~and(abs(nodes(:,:,1)-ctr) <= R, abs(nodes(:,:,2)-ctr) <= R);
        elseif dim == 3
            isFree = ~and(abs(nodes(:,:,:,1)-ctr) <= R, and(abs(nodes(:,:,:,2)-ctr) <= R,abs(nodes(:,:,:,3)-ctr) <= R));
        end
    elseif strcmpi(geometry,'squareSlowdown')
        mDim = m*ones(1,dim);
        isFree = true(mDim);
    else
        error('Geometry must be circle or square.');
    end
end

isFree = double(isFree);
numFree = sum(isFree(:) > 0);
isFree(isFree > 0) = 1:numFree;

nodes = reshape(nodes,m^dim,dim);
nodes = nodes(isFree(:) > 0,:);

time = toc;
fprintf('Calculated available sites in %.1f seconds.\n',time);

end

