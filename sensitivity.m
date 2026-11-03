% SENSITIVITY

model = sg_model_init();
outputDir = fileparts(mfilename('fullpath'));
lambdaGrid = linspace(250,7500,70);
qs = linspace(0.03,1,300);
qStar = zeros(numel(lambdaGrid),2);
maxASE = zeros(numel(lambdaGrid),2);
for i = 1:numel(lambdaGrid)
    model.lambda = lambdaGrid(i);
    for j = 1:2
        Ns = [0,128];
        cv = sg_coverage(qs,0.05,model,Ns(j),1200,true);
        ase = 2*model.lambda*qs.*cv;
        [maxASE(i,j), idx] = max(ase);
        qStar(i,j) = qs(idx);
    end
end
fprintf('Pair-density sensitivity: [lambda, q* no RIS, ASE no RIS, q* RIS, ASE RIS]\n');
for ix = [3,22,46,70]   
    fprintf('%8.1f  %.5f  %.5f  %.5f  %.5f\n', ...
            lambdaGrid(ix),qStar(ix,1),maxASE(ix,1), ...
            qStar(ix,2),maxASE(ix,2));
end
fig = figure('Color','w','Position',[100 100 620 360]);
plot(lambdaGrid,qStar(:,2),'-','LineWidth',3, ...
     'DisplayName','FD + RIS (128 elements)');
hold on;
plot(lambdaGrid,qStar(:,1),'--','LineWidth',3, ...
     'DisplayName','FD without RIS');
xlabel('Potential D2D pair density (pairs/km^2)');
ylabel('ASE-optimal channel access probability');
ylim([0,1.04]); grid on; box on;
legend('show','Location','best');
sg_save_figure(fig,'fig_density',outputDir,160);
