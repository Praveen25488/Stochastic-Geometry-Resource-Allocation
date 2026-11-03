function m = sg_mean_mark(d, s, model)
%SG_MEAN_MARK Exponential-fading + LOS/NLOS + antenna-gain mark average.
%  D can be a column vector or matrix. The result has the same size.

pLOS = exp(-model.kappa*d);
Llos = (1+d).^(-model.alphaL);
Lnlos = (1+d).^(-model.alphaN);
m = zeros(size(d));
for j = 1:numel(model.GG)
    g = model.GG(j);
    m = m + model.probs(j)*(pLOS./(1+s*g*Llos) + ...
                           (1-pLOS)./(1+s*g*Lnlos));
end
end
