% VALIDATE_MC.M
% Spatial Monte Carlo verification.
% The conditional fading, blockage and antenna marks are averaged
% analytically with SG_MEAN_MARK; PPP locations are Monte Carlo sampled.


model = sg_model_init();
rng(20261014, 'twister');
Rcut = 700;
nrun = 2000;
lambdaM2 = model.lambda/1e6;

for q = [0.4, 0.8]
    samples = zeros(nrun,1);
    [g,w] = sg_link_states(model,128,model.lambdaR);
    for i = 1:nrun
        meanN = lambdaM2*q*pi*Rcut^2;
        n = sg_poisson_rand(meanN);
        rr = Rcut*sqrt(rand(n,1));
        ang = 2*pi*rand(n,1);
        angleEnd = 2*pi*rand(n,1);
        x = rr.*cos(ang);
        y = rr.*sin(ang);
        rr2 = hypot(x+model.r0*cos(angleEnd), ...
                    y+model.r0*sin(angleEnd));
        samplePC = 0;
        for j = 1:numel(g)
            if w(j)==0, continue; end
            s = model.T/g(j);
            m1 = sg_mean_mark(rr,s,model);
            m2 = sg_mean_mark(rr2,s,model);
            % log of a product avoids numerical underflow.
            logProd = sum(log(m1)+log(m2));
            samplePC = samplePC + w(j)* ...
                exp(-model.T*(model.beta+model.sigma/model.P)/g(j) + logProd);
        end
        samples(i) = samplePC;
    end
    theory = sg_coverage(q,model.P,model,128,model.lambdaR,true);
    meanMC = mean(samples);
    standardError = std(samples,0)/sqrt(nrun);
    ci = meanMC + [-1,1]*1.96*standardError;
    fprintf('q=%.1f: analysis=%.5f; MC=%.5f; 95%% CI=[%.5f, %.5f]\n', ...
        q,theory,meanMC,ci(1),ci(2));
end
