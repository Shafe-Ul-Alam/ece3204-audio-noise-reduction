function [a, d] = dwt_periodic(x, h)
% DWT_PERIODIC  One level of the orthogonal discrete wavelet transform.
%
%   [a, d] = dwt_periodic(x, h) splits x (even length) into
%     a : approximation (low-pass) coefficients, length N/2
%     d : detail (high-pass) coefficients,       length N/2
%   using the low-pass scaling filter h and the quadrature-mirror high-pass
%   filter g[n] = (-1)^n h[L-1-n], with periodic signal extension.
%
%   a[k] = sum_n h[n] x[(2k+n) mod N]      (filter + downsample by 2)
%   d[k] = sum_n g[n] x[(2k+n) mod N]
%
%   Implemented from scratch so the project does not need the Wavelet Toolbox.

x = x(:); h = h(:);
N = length(x); L = length(h);
g = qmf(h);

k   = (0:N/2-1)';
idx = mod(2*k + (0:L-1), N) + 1;     % (N/2 x L) matrix of sample indices
X   = x(idx);
a   = X * h;
d   = X * g;
end

function g = qmf(h)
L = length(h);
g = flipud(h) .* (-1).^(0:L-1)';
end
