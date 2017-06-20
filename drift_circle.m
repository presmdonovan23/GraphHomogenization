function val = drift_circle(x,j,rho)
% only makes sense for circular obstuctions

ctr = .5;
K1 = 25;
K2 = 10;
R = rho/2;
dist2ob = @(x) sqrt(sum((x-ctr).^2,2))-R;

xNormalized = (x(:,j) - ctr)./sqrt(sum((x-ctr).^2,2));

val = K1*xNormalized.*exp(-K2*dist2ob(x));
val(~isfinite(val)) = 0;

end