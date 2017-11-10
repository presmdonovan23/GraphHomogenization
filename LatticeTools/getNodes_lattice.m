function [nodes, isFree] = getNodes_lattice( latticeGeo, verbose )
% assumes nodes live on lattice
% assumes cell has length 1

if nargin < 2 || isempty(verbose)
    verbose = 1;
end

tic

dim = latticeGeo.dim;
m = latticeGeo.m;
h = latticeGeo.h;
obRad = latticeGeo.obRad;
name = latticeGeo.name;
obCtr = latticeGeo.obCtr;
specialSetting = latticeGeo.specialSetting;

if dim == 2
    obCtr3(1,1,:) = obCtr;
elseif dim == 3
    obCtr3(1,1,1,:) = obCtr;
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
nodes = (nodes - .5)*h;

if obRad == 0 || strcmpi(specialSetting,'m2_slowOneSite')
    mDim = m*ones(1,dim);
    isFree = true(mDim);
else
    if strcmpi(name,'circle')
        dist2ctr2 = sum((nodes - obCtr3).^2,dim+1);
        dist2ctr2 = round(dist2ctr2,m); % eliminate numerical error that could make obstructed sites non-symmetrical
        isFree = dist2ctr2 > obRad^2;
    elseif strcmpi(name,'square')
        % ** could be some numerical error for small m where sites are on boundary of obstructed region
        isFree = ~all(abs(nodes - obCtr3) <= obRad,dim+1);
    end
end

isFree = double(isFree);
numFree = sum(isFree(:) > 0);
isFree(isFree > 0) = 1:numFree;

nodes = reshape(nodes,m^dim,dim);
nodes = nodes(isFree(:) > 0,:);

time = toc;
if verbose
    fprintf('\tCalculated available sites in %.1f seconds.\n',time);
end

end