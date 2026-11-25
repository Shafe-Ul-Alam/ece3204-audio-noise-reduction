function m = compute_metrics(clean, noisy, y, fs, cfg)
% COMPUTE_METRICS  Objective quality measures of a processed signal y.
%
%   m.snr_in   : SNR of the noisy input (dB)
%   m.snr_out  : SNR of the processed output (dB)
%   m.dsnr     : SNR improvement = snr_out - snr_in (dB)
%   m.segsnr   : segmental SNR of the output (dB) - closer to perceived quality
%   m.lsd      : log-spectral distance to the clean signal (dB) - lower is better,
%                penalises spectral distortion / removed speech components
%   m.corr     : correlation coefficient between clean and output

clean = clean(:); noisy = noisy(:); y = y(:);

m.snr_in  = snr_db(clean, noisy);
m.snr_out = snr_db(clean, y);
m.dsnr    = m.snr_out - m.snr_in;
m.segsnr  = seg_snr(clean, y, cfg.seg_len);
m.lsd     = log_spectral_distance(clean, y, cfg.seg_len);
c = corrcoef(clean, y);
m.corr    = c(1, 2);
end

% --------------------------------------------------------------------------
function s = snr_db(ref, est)
s = 10*log10(sum(ref.^2) / sum((ref - est).^2));
end

function s = seg_snr(ref, est, L)
% Segmental SNR: average of frame SNRs, each clamped to [-10, 35] dB.
% Silent frames (energy < 40 dB below the loudest frame) are skipped.
nfr = floor(length(ref) / L);
R = reshape(ref(1:nfr*L), L, nfr);
E = reshape(est(1:nfr*L), L, nfr);
Pr = sum(R.^2, 1);
Pe = sum((R - E).^2, 1) + eps;
active = Pr > max(Pr) * 1e-4;
fsnr = 10*log10(Pr(active) ./ Pe(active) + eps);
fsnr = min(max(fsnr, -10), 35);
s = mean(fsnr);
end

function d = log_spectral_distance(ref, est, L)
% RMS difference of log power spectra (dB), averaged over active frames,
% with a 50 dB dynamic-range floor per frame.
nfr = floor(length(ref) / L);
w = 0.54 - 0.46*cos(2*pi*(0:L-1)'/(L-1));          % Hamming window
R = reshape(ref(1:nfr*L), L, nfr) .* w;
E = reshape(est(1:nfr*L), L, nfr) .* w;
PR = abs(fft(R)).^2; PR = PR(1:L/2+1, :);
PE = abs(fft(E)).^2; PE = PE(1:L/2+1, :);
active = sum(PR, 1) > max(sum(PR, 1)) * 1e-4;
PR = PR(:, active); PE = PE(:, active);
% limit the dynamic range to 50 dB below each frame's peak so that bins
% with (almost) no speech energy do not dominate the distance
floor_ = max(PR, [], 1) * 1e-5;
D = 10*log10(max(PR, floor_)) - 10*log10(max(PE, floor_));
d = mean(sqrt(mean(D.^2, 1)));
end
