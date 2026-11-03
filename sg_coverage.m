function pc = sg_coverage(q, p, model, N, lambdaR, isFD)
%SG_COVERAGE SINR coverage probability from paired-FD PPP quadrature.

if nargin < 4, N = 128; end
if nargin < 5, lambdaR = model.lambdaR; end
if nargin < 6, isFD = true; end
if any(p(:) <= 0) || any(q(:) < 0) || any(q(:) > 1)
    error('sg_coverage:Range', 'Require p > 0 and 0 <= q <= 1.');
end
if ~isscalar(q) && ~isscalar(p) && ~isequal(size(q),size(p))
    error('sg_coverage:Shape', 'Non-scalar q and p must have matching sizes.');
end
pc = zeros(size(q+p));
[g, w] = sg_link_states(model, N, lambdaR);
if isFD
    kvals = model.kPair;
    b = model.beta;
else
    kvals = model.kSingle;
    b = 0;
end
for t = 1:numel(g)
    if w(t) == 0, continue; end
    [~, ix] = min(abs(model.kGains - g(t)));
    if abs(model.kGains(ix)-g(t)) > 1e-13
        error('sg_coverage:MissingGain', ...
              'RIS gain N=%g is outside the precomputed table.', N);
    end
    pc = pc + w(t).*exp(-model.T*(b+model.sigma./p)/g(t) ...
                   -(model.lambda/1e6).*q.*kvals(ix));
end
end
