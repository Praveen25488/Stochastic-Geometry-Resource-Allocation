% SI_SENSITIVITY

model = sg_model_init();
outputDir = fileparts(mfilename('fullpath'));
qq = linspace(0.03,1,160);
betaDb = linspace(-85,-35,91);
maxASE = zeros(numel(betaDb),3);
qStar = zeros(numel(betaDb),3);
Ns = [128,0,0];
FDcases = [true,true,false];
for i = 1:numel(betaDb)
    model.beta = 10^(betaDb(i)/10);
    for j = 1:3
        fd = FDcases(j);
        cv = sg_coverage(qq,0.05,model,Ns(j),model.lambdaR,fd);
        ase = (1+double(fd))*model.lambda*qq.*cv;
        [maxASE(i,j), idx] = max(ase);
        qStar(i,j) = qq(idx);
    end
end
fprintf('Residual-SI results: beta(dB), ASE [FD+RIS FD HD], q* [FD+RIS FD HD]\n');
for target = [-80,-60,-50,-40,-35]
    [~, k] = min(abs(betaDb-target));
    fprintf('%4d   ASE %9.3f %9.3f %9.3f   q* %.3f %.3f %.3f\n', ...
        target,maxASE(k,1),maxASE(k,2),maxASE(k,3), ...
        qStar(k,1),qStar(k,2),qStar(k,3));
end
fig = figure('Color','w','Position',[100 100 620 360]);
plot(-betaDb,maxASE(:,1),'-','LineWidth',3, ...
     'DisplayName','FD + RIS (128 elements)'); hold on;
plot(-betaDb,maxASE(:,2),'--','LineWidth',3, ...
     'DisplayName','FD without RIS');
plot(-betaDb,maxASE(:,3),':','LineWidth',3, ...
     'DisplayName','HD without RIS');
xlabel('Residual self-interference suppression (dB)');
ylabel('Peak ASE (bit/s/Hz/km^2)');
legend('show','Location','best'); grid on; box on;
sg_save_figure(fig,'fig_si',outputDir,160);
