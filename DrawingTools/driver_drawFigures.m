saveOn = 0;

%% setting
rho = .45;
m = 8;
geometry = 'circle';
drawObs = 1;
ctr = .5;
drawSetting(rho,m,ctr,geometry,drawObs,saveOn)

%% path length effects
rho = .5;
geometry = 'square';
mVals = [4 8 16 32];
drawObs = 1;

drawPathLengthEffects(rho,geometry,mVals,drawObs,saveOn)
%% homog theory lim
rho = .45;
geometry = 'circle';
m = 16;
epsInvVals = [1 2 4 8];
drawObs = 0;

drawHomogLimit(rho,geometry,m,epsInvVals,drawObs,saveOn)
%% drift field
rho = .45;
m = 8;
geometry = 'circle';
alpha = 1;
drawObs = 1;

drawDriftField(rho,m,geometry,alpha,drawObs,saveOn);
 
