load('compTimeTests.mat')

mVals = [ghParams.m];
n = length(mVals);

% qr decomposition

clear time_qr res_qr
for i = 1:8
    i
    L = ghInput(i).L;
    
    tic;
    [Q,~,~] = qr(L);
    pi0 = Q(:,end);
    pi0 = pi0./sum(pi0);
    
    time_qr(i) = toc;
    res_qr(i) = norm(L'*pi0);
    
end
save('compTimeResults.mat','time_*','res_*');
% svd decomposition

clear time_svd res_svd
for i = 1:7
    i
    L = ghInput(i).L;
    
    tic;
    
    [~,~,V] = svds(L',1,'smallest');
    pi0 = V;
    pi0 = pi0./sum(pi0);
    
    time_svd(i) = toc;
    res_svd(i) = norm(L'*pi0);
    
end
save('compTimeResults.mat','time_*','res_*');
% eigs

clear time_eigs res_eigs
for i = 1:n
    i
    L = ghInput(i).L;
    
    tic;
    
    [pi0,~] = eigs(L',1,0);
    pi0 = pi0./sum(pi0);
    
    time_eigs(i) = toc;
    res_eigs(i) = norm(L'*pi0);
    
end
save('compTimeResults.mat','time_*','res_*');
% inverse iteration

clear time_ii res_ii
for i = 1:n
    i
    L = ghInput(i).L;
    
    tic;
    
    [LL,U,P,Q,R] = lu(L');
    nNodes = size(L,1);
    pi0 = randn(nNodes,1);
    pi0 = pi0./max(abs(pi0));
    
    pi0 = Q * (U \ (LL \ (P * (R \ pi0))));
    pi0 = pi0./sum(pi0);
    pi0 = Q * (U \ (LL \ (P * (R \ pi0))));
    pi0 = pi0./sum(pi0);

    time_ii(i) = toc;
    res_ii(i) = norm(L'*pi0);
    
end
save('compTimeResults.mat','time_*','res_*');
clearvars -except time_* res_* mVals

%% plot
load('compTimeTests.mat')
load('compTimeResults.mat')

ghParams(1) = [];
ghInput(1) = [];
time_qr(1) = [];
res_qr(1) = [];
time_svd(1) = [];
res_svd(1) = [];
time_eigs(1) = [];
res_eigs(1) = [];
time_ii(1) = [];
res_ii(1) = [];

mVals = [ghParams.m];
n = length(mVals);

fh = figure
hold on
%plot(mVals(1:length(time_qr)),time_qr);
%plot(mVals(1:length(time_svd)),time_svd);
%plot(mVals(1:length(time_eigs)),time_eigs);
%plot(mVals(1:length(time_ii)),time_ii);
log2m = log2(mVals);
plot(log2m(1:length(time_qr)),log10(time_qr),'linewidth',3);
plot(log2m(1:length(time_svd)),log10(time_svd),'linewidth',3);
plot(log2m(1:length(time_eigs)),log10(time_eigs),'linewidth',3);
plot(log2m(1:length(time_ii)),log10(time_ii),'linewidth',3);
legend('QR Decomposition','SVD','\texttt{eigs}','Inverse Iteration','location','southeast')
ca = gca;
ca.FontSize = 18;
ca.XTick = [log2(mVals)];
for i = 1:length(ca.XTickLabel)
    ca.XTickLabel{i} = sprintf('2^{-%d}',log2(mVals(i)));
end

%yticks = [1e-5 1e-4 1e-3 1e-2 1e-1 1 10 20];
%ca.YTick = log10([1e-5 1e-4 1e-3 1e-2 1e-1 1 10 20]);
for i = 1:length(ca.YTick)
    expStr = str2double(ca.YTickLabel{i});
    
    if expStr == 0
        ca.YTickLabel{i} = '1';
    elseif expStr == 1
        ca.YTickLabel{i} = sprintf('10');
    else
        ca.YTickLabel{i} = sprintf('10^{%d}',expStr);
    end
end
xlabel('Path length');
ylabel('Computation time (sec)')

dirname = '/Users/prestondonovan/Documents/School/Research/Thesis/includes/chapter2/';

mySaveFig([dirname 'compTimeStatDist'],fh,'fig')
mySaveFig([dirname 'compTimeStatDist'],fh,'png')
mySaveFig([dirname 'compTimeStatDist'],fh,'eps')