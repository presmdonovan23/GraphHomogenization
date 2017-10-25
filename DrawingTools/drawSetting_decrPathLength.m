%% Plot graph setting w/ square in upper right corner (m = 2,4,8,16)
saveOn = 0;

fh = figure;
i = 1;
ctr = .75;
geometryName = 'square';
drawObs = 1;
drawEdges = 1;
drawRates = 0;
loc = [];

for m = [2 4 8 16]
    sph(i) = subplot(2,2,i);
    if i == 2 || i == 4
        sph(i).Position = sph(i).Position - [.15 .1 -.15 -.15];
    else
        sph(i).Position = sph(i).Position - [.1 .1 -.15 -.15];
    end
    rho = .5;
    %m = 2;
    
    drawCell_lattice(rho,m,ctr,geometryName,drawObs,drawEdges,drawRates,saveOn,loc,fh);
    i=i+1;
end

if saveOn
    dirname = '/Users/prestondonovan/Documents/School/Research/My Notes & Papers/[Donovan,Rathinam]_Graph_Homogenization Publication/publication/figures/';

    mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'fig')
    mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'png')
    mySaveFig([dirname 'pathLengthEffectsGraph'],fh,'eps')
end