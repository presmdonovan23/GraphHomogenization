function val = rate_square_stick(x,h,D0,s)
% s = side length of square

d = .75*h; % will only affect sites at boundary
%d = 1/8;
penalty = .95;

diffusionRate = D0/h^2;

dist2ob_1norm = abs(x - .5) - s/2;
nearBoundary = prod(dist2ob_1norm < d,2);

% if at boundary, reduce rate by 95%
val = diffusionRate - nearBoundary*penalty*diffusionRate;

end