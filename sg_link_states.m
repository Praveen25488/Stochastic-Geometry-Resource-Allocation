function [g, w] = sg_link_states(model, N, lambdaR)
%SG_LINK_STATES Path-selection outcomes and association weights.

%   Column arrays G and W represent the four LOS/NLOS/RIS states.

if nargin < 3, lambdaR = model.lambdaR; end
[gL, gN, gR] = sg_gain(model, N);
if N > 0
    avail = 1-exp(-(lambdaR/1e6)*pi*model.R^2*exp(-model.kappa*model.r0));
else
    avail = 0;
end
pL = exp(-model.kappa*model.r0);
g = [gL; max(gL,gR); max(gN,gR); gN];
w = [pL*(1-avail); pL*avail; (1-pL)*avail; (1-pL)*(1-avail)];
end
