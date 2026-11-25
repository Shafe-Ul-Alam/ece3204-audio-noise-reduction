function [P, f] = psd_periodogram(x, fs, nfft)
% PSD_PERIODOGRAM  One-sided periodogram PSD estimate (rectangular window).
%
%   P(f) = |X(f)|^2 / (fs * N), doubled for 0 < f < fs/2  [units: power/Hz]
%   The whole signal is one segment -> fine frequency resolution but very high
%   variance (the estimate does not get smoother as N grows).

x = x(:);
N = length(x);
if nargin < 3 || isempty(nfft), nfft = 2^nextpow2(N); end
nfft = max(nfft, N);

X = fft(x, nfft);
P = abs(X).^2 / (fs * N);
P = P(1:floor(nfft/2)+1);
P(2:end-1) = 2 * P(2:end-1);
f = (0:floor(nfft/2))' * fs / nfft;
end
