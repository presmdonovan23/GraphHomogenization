%% Continuous homogenization vs discrete homogenization
saveOn = 0;

load('Results_2d_square/results_2017_08_05_10_35_27.mat','results_homog')

x_homog = [1:9];
y_homog = [2/3 results_homog.Deff];
x_cont = linspace(1,9,100);
y_cont = .77*ones(1,100);

fh = figure;

hold on
plot(x_cont,y_cont,'r-','linewidth',2)
plot(x_homog,y_homog,'b.-','markersize',20,'linewidth',2)

axis([1 9 .66 .78]);
axis square

xlabel('Path Length (h)');
ylabel('D_e');

lh = legend('Continuous Homogenization','Graph Homogenization','location','southeast');

ca = gca;
ca.FontSize = 18;
ca.XTick = [1:9];
for i = 1:9
    ca.XTickLabel{i} = sprintf('2^{-%d}',i);
end

if saveOn
    dirname = '';

    mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'fig')
    mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'png')
    mySaveFig([dirname 'pathLengthEffectsDeff'],fh,'eps')
end