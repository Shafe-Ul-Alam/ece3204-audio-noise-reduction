function [X, f] = fft_db(x, fs)
% FFT_DB  One-sided FFT magnitude spectrum in dB (normalised to N).
x = x(:);
N = length(x);
nfft = 2^nextpow2(N);
Xc = fft(x, nfft) / N;
X = 20*log10(abs(Xc(1:nfft/2+1)) + 1e-12);
f = (0:nfft/2)' * fs / nfft;
end
