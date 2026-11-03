function k = sg_poisson_rand(mu)
%SG_POISSON_RAND Poisson random count without Statistics Toolbox.
%   For mu>=30 uses the transformed-rejection method (PTRS).
%   For smaller means uses the classic product-of-uniforms algorithm.

if ~isscalar(mu) || mu<0 || ~isfinite(mu)
    error('sg_poisson_rand:InvalidMu','mu must be a finite nonnegative scalar.');
end
if mu == 0, k = 0; return; end
if mu < 30
    lim = exp(-mu);
    prodVal = 1;
    k = 0;
    while true
        k = k+1;
        prodVal = prodVal*rand();
        if prodVal <= lim
            k = k-1;
            return;
        end
    end
end
b = 0.931 + 2.53*sqrt(mu);
a = -0.059 + 0.02483*b;
invAlpha = 1.1239 + 1.1328/(b-3.4);
vr = 0.9277 - 3.6224/(b-2);
while true
    u = rand()-0.5;
    v = rand();
    us = 0.5-abs(u);
    if us <= 0, continue; end
    k = floor((2*a/us + b)*u + mu + 0.43);
    if k >= 0 && us >= 0.07 && v <= vr
        return;
    end
    if k < 0 || (us < 0.013 && v > us)
        continue;
    end
    lhs = log(v) + log(invAlpha) - log(a/(us^2)+b);
    rhs = -mu + k*log(mu) - gammaln(k+1);
    if lhs <= rhs
        return;
    end
end
end
