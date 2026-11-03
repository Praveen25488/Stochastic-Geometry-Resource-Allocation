function model = sg_model_init()

persistent baseModel
if ~isempty(baseModel)
    model = baseModel;
    return;
end

model.r0 = 40;
model.kappa = 0.015;
model.alphaL = 2.2;
model.alphaN = 4.0;
model.side = 0.05;
model.theta = pi/6;
model.prob0 = model.theta/(2*pi);
model.GG = [1, model.side, model.side^2];
p0 = model.prob0;
model.probs = [p0^2, 2*p0*(1-p0), (1-p0)^2];

model.T = 1;
model.beta = 1e-6;
model.sigma = 1e-8;
model.P = 0.05;
model.lambda = 2500;             % D2D pairs / km^2
model.lambdaR = 1200;            % RIS panels / km^2
model.R = 20;                    % RIS availability radius (m)
model.etaR = 0.003;

r = [0; logspace(-4, log10(3000), 2000)'];
u = (0:71)*(2*pi/72);        % endpoint angular quadrature
r_other = sqrt(r.^2 + model.r0^2 + 2*(r*cos(u))*model.r0);

[gL, gN] = sg_gain(model, 0);
[~, ~, gR64] = sg_gain(model, 64);
[~, ~, gR128] = sg_gain(model, 128);
[~, ~, gR256] = sg_gain(model, 256);
gains = unique([gL, gN, gR64, gR128, gR256], 'stable');

pairs = zeros(size(gains));
singles = zeros(size(gains));
for ix = 1:numel(gains)
    s = model.T / gains(ix);
    m1 = sg_mean_mark(r, s, model);
    m2 = mean(sg_mean_mark(r_other, s, model), 2);
    pairs(ix) = 2*pi*sg_simpson(r, r.*(1-m1.*m2));
    singles(ix) = 2*pi*sg_simpson(r, r.*(1-m1));
end

model.kGains = gains;
model.kPair = pairs;
model.kSingle = singles;
baseModel = model;
end
