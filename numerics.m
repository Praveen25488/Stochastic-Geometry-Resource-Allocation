model = sg_model_init();
outputDir = fileparts(mfilename('fullpath'));

% Fig. 1: area spectral efficiency versus channel access probability.
q = linspace(0.01, 1, 160);
fig = figure('Color','w','Position',[100 100 620 380]);
hold on;
Ncases = [128, 0, 0];
FDcases = [true, true, false];
styles = {'-', '--', ':'};
labels = {'FD + RIS (128 elements)', 'FD without RIS', 'HD without RIS'};
for ix = 1:3
    isFD = FDcases(ix);
    pc = sg_coverage(q, model.P, model, Ncases(ix), model.lambdaR, isFD);
    factor = 1 + double(isFD);
    ase = factor*model.lambda*q.*pc;
    [val, optix] = max(ase);
    fprintf('%s: max ASE %.3f, q*=%.3f, ASE q=1 %.3f, coverage*=%.4f\n', ...
            labels{ix}, val, q(optix), ase(end), pc(optix));
    plot(q, ase, styles{ix}, 'LineWidth',3, 'DisplayName',labels{ix});
end
xlabel('Channel access probability, q');
ylabel('ASE (bit/s/Hz/km^2)');
legend('show','Location','best'); grid on; box on;
sg_save_figure(fig, 'fig_ase', outputDir, 180);

% Fig. 2: coverage versus RIS density.
fig = figure('Color','w','Position',[100 100 620 380]);
hold on;
lamr = linspace(0,3000,85);
Ns = [64,128,256];
styles = {'--','-',':'};
for ix = 1:numel(Ns)
    N = Ns(ix);
    vals = arrayfun(@(lr) sg_coverage(0.6,model.P,model,N,lr,true),lamr);
    plot(lamr,vals,styles{ix},'LineWidth',3, ...
         'DisplayName',sprintf('N = %d',N));
end
fprintf('RIS compare N=128 q=0.6:\n');
for lr = [0,500,1200,3000]
    fprintf('  lambdaR=%4d: coverage %.5f\n',lr, ...
            sg_coverage(0.6,model.P,model,128,lr,true));
end
xlabel('RIS density (surfaces/km^2)');
ylabel('SINR coverage probability at T = 0 dB');
legend('show','Location','best'); grid on; box on;
sg_save_figure(fig, 'fig_ris', outputDir, 180);

% EE maximisation on the same finite q-p grid as the original script.
qs = linspace(0.05,1,101);
ps = logspace(log10(0.0001),log10(0.2),101);
Pc = 0.05; xi = 0.35; PR = 0.01;
bestEE = -inf; bestQ = NaN; bestP = NaN;
for pw = ps
    cv = sg_coverage(qs,pw,model,128,model.lambdaR,true);
    ase = 2*model.lambda*qs.*cv;
    denom = 2*model.lambda*qs.*(pw/xi+Pc)+model.lambdaR*PR;
    ee = ase./denom;
    [optEE, optix] = max(ee);
    if optEE > bestEE
        bestEE = optEE;
        bestQ = qs(optix);
        bestP = pw;
    end
end
fprintf('EE grid optimum: EE=%.8f, q*=%.6f, p*=%.10f W\n', ...
        bestEE,bestQ,bestP);
for pw = [1e-3,1e-2,0.05,0.1]
    cv = sg_coverage(0.6,pw,model);
    fprintf('q=0.6, p=%g: coverage %.9f, ASE %.5f\n', ...
            pw,cv,2*model.lambda*0.6*cv);
end
[g,w] = sg_link_states(model,128,model.lambdaR);
fprintf('Link-state gain and weight [gain, probability]:\n');
disp([g,w]);
fprintf('K-pair, K-single gain table:\n');
disp([model.kGains(:),model.kPair(:),model.kSingle(:)]);
