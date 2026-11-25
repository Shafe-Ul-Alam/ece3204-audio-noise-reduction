function [y, w_hist] = nlms_hum_canceller(x, fs, cfg)
% NLMS_HUM_CANCELLER  Adaptive noise canceller for mains hum (adaptive notch).
%
%   Primary input  d[n] = speech + hum
%   Reference r[n]      = [cos(w1 n) sin(w1 n) ... cos(wK n) sin(wK n)],
%                         wk = 2*pi*k*f0/fs  (f0 = 50 Hz and its harmonics)
%   The NLMS filter learns the amplitude/phase of each hum harmonic:
%       yhat[n] = w' r[n]                   (hum estimate)
%       e[n]    = d[n] - yhat[n]            (cleaned output)
%       w      <- w + mu * e[n] r[n] / (eps + r' r)
%   The filter tracks slow changes in hum amplitude/phase automatically,
%   which is what makes it "adaptive". It does nothing against white noise
%   or hiss - those are not correlated with the reference.

x = x(:);
N = length(x);
K = cfg.nlms_harm;
n = (0:N-1)';

R = zeros(N, 2*K);
for k = 1:K
    R(:, 2*k-1) = cos(2*pi*k*cfg.hum_f0*n/fs);
    R(:, 2*k)   = sin(2*pi*k*cfg.hum_f0*n/fs);
end

w  = zeros(2*K, 1);
y  = zeros(N, 1);
mu = cfg.nlms_mu;
store = nargout > 1;
if store, w_hist = zeros(N, 2*K); end

for i = 1:N
    r    = R(i, :)';
    e    = x(i) - w' * r;
    w    = w + (mu / (1e-6 + r' * r)) * e * r;
    y(i) = e;
    if store, w_hist(i, :) = w'; end
end
end
