function x = idwt_periodic(a, d, h)
% IDWT_PERIODIC  Inverse of DWT_PERIODIC (perfect reconstruction).
%
%   x[m] = sum_k ( a[k] h[m-2k] + d[k] g[m-2k] ),  indices taken mod N
%   (upsample by 2 + synthesis filtering). Because the filter bank is
%   orthogonal, the synthesis filters are the time-reverses of the analysis
%   filters, which is exactly what this transpose operation implements.

a = a(:); d = d(:); h = h(:);
L = length(h); N = 2*length(a);
g = flipud(h) .* (-1).^(0:L-1)';

k   = (0:N/2-1)';
idx = mod(2*k + (0:L-1), N) + 1;            % same index map as the analysis
contrib = a * h.' + d * g.';                % (N/2 x L) contributions
x = accumarray(idx(:), contrib(:), [N 1]);
end
