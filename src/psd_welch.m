function [P, f] = psd_welch(x, fs, win_len, overlap, nfft)
% PSD_WELCH  Welch's averaged, modified periodogram (Hamming window).
%
%   The signal is split into overlapping segments of win_len samples, each is
%   windowed, a periodogram is computed and all periodograms are averaged.
%   Averaging K segments reduces the variance by roughly a factor of K, at the
%   cost of frequency resolution (~ fs / win_len).

x = x(:);
if nargin < 5 || isempty(nfft), nfft = win_len; end

w    = 0.54 - 0.46*cos(2*pi*(0:win_len-1)'/(win_len-1));   % Hamming
U    = sum(w.^2);                                            % window power
hop  = round(win_len * (1 - overlap));
nseg = floor((length(x) - win_len) / hop) + 1;

P = zeros(nfft, 1);
for k = 1:nseg
    seg = x((k-1)*hop + (1:win_len)) .* w;
    P   = P + abs(fft(seg, nfft)).^2;
end
P = P / (nseg * fs * U);

P = P(1:floor(nfft/2)+1);
P(2:end-1) = 2 * P(2:end-1);
f = (0:floor(nfft/2))' * fs / nfft;
end
