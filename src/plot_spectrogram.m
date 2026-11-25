function plot_spectrogram(x, fs, ttl)
% PLOT_SPECTROGRAM  Short-time Fourier transform image (no toolbox needed).
%   Hamming window of 512 samples (32 ms at 16 kHz), 75 % overlap.
x = x(:);
L = 512; hop = 128; nfft = 512;
w = 0.54 - 0.46*cos(2*pi*(0:L-1)'/(L-1));
nfr = floor((length(x) - L) / hop) + 1;
S = zeros(nfft/2+1, nfr);
for k = 1:nfr
    seg = x((k-1)*hop + (1:L)) .* w;
    X = fft(seg, nfft);
    S(:, k) = 20*log10(abs(X(1:nfft/2+1)) + 1e-10);
end
t = ((0:nfr-1)*hop + L/2) / fs;
f = (0:nfft/2) * fs / nfft;
imagesc(t, f/1000, S);
axis xy;
mx = max(S(:));
caxis([mx-80 mx]);
colormap(jet);
colorbar;
xlabel('Time (s)'); ylabel('Frequency (kHz)');
title(ttl);
end
