function I = sg_simpson(x, y)
%SG_SIMPSON Composite Simpson integration on NON-UNIFORM radial grid.

x = x(:);
y = y(:);
if numel(x) ~= numel(y) || numel(x) < 3
    error('sg_simpson:InvalidInput', 'X and Y must have >= 3 elements.');
end
n = numel(x);
if mod(n,2) ~= 1
    error('sg_simpson:EvenN', 'Use an odd number of quadrature nodes.');
end
a = 1:2:(n-2); b = a+1; c = a+2;
h0 = x(b)-x(a);
h1 = x(c)-x(b);
h = h0+h1;
term = (h/6).*( (2-h1./h0).*y(a) + ...
                 (h.^2./(h0.*h1)).*y(b) + ...
                 (2-h0./h1).*y(c) );
I = sum(term);
end
