function [P, f, a, E] = ar_psd_estimate(x, p, nfft, fs, method)
% AR_PSD_ESTIMATE  Parametric (autoregressive) PSD estimate, implemented from scratch.
%
%   Model: x(n) = -a1 x(n-1) - ... - ap x(n-p) + e(n),  var(e) = E
%   PSD  : P(f) = E / ( fs |1 + sum_k a_k e^{-j 2 pi f k / fs}|^2 ), one-sided
%
%   method 'yule' : Yule-Walker - solve the Toeplitz autocorrelation equations
%                   R a = -r  (same model as MATLAB pyulear)
%   method 'burg' : Burg - choose each reflection coefficient to minimise the
%                   forward + backward prediction error (same as pburg);
%                   better resolution for short records, always stable.

x = x(:) - mean(x);
N = length(x);
switch lower(method)
    case 'yule'
        r = zeros(p+1, 1);
        for k = 0:p
            r(k+1) = sum(x(1:N-k) .* x(1+k:N)) / N;     % biased autocorrelation
        end
        R = toeplitz(r(1:p));
        a = [1; -(R \ r(2:p+1))];
        E = r(1) + sum(a(2:end) .* r(2:p+1));
    case 'burg'
        ef = x; eb = x; a = 1; E = (x' * x) / N;
        for m = 1:p
            efp = ef(2:end); ebp = eb(1:end-1);
            k   = -2 * (ebp' * efp) / (efp' * efp + ebp' * ebp);
            ef  = efp + k * ebp;
            eb  = ebp + k * efp;
            a   = [a; 0] + k * [0; flipud(a)];
            E   = (1 - k^2) * E;
        end
    otherwise
        error('method must be ''yule'' or ''burg''');
end

A = fft(a, nfft);
P = E ./ (fs * abs(A).^2);
P = P(1:floor(nfft/2)+1);
P(2:end-1) = 2 * P(2:end-1);
f = (0:floor(nfft/2))' * fs / nfft;
end
