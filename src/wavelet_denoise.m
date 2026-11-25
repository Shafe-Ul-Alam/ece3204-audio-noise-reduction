function [y, info] = wavelet_denoise(x, cfg)
% WAVELET_DENOISE  DWT-based denoising (wavelet chosen by cfg.wavelet, default db4).
%
%   Steps
%     1. Multi-level DWT (cfg.wav_levels levels) -> a_J, d_J, ..., d_1
%     2. Estimate the noise level sigma with the robust median estimator
%           sigma = median(|d_1|) / 0.6745   (finest level = mostly noise)
%     3. Threshold per level, three rules (cfg.wav_rule):
%        'bayes'  : BayesShrink  T_j = sigma^2 / sigma_x,j,
%                   sigma_x,j = sqrt(max(var(d_j) - sigma^2, 0))
%                   -> adapts to the signal energy in every band (default)
%        'global' : VisuShrink universal threshold T = sigma*sqrt(2 ln N)
%        'level'  : universal threshold with a separate sigma per level
%     4. Soft (or hard) thresholding of the detail coefficients only
%     5. Inverse DWT
%
%   info returns coefficients before/after thresholding (for plots).

x = x(:);
N0 = length(x);
J  = cfg.wav_levels;

% pad to a multiple of 2^J (symmetric extension avoids an artificial jump)
Np = ceil(N0 / 2^J) * 2^J;
xp = [x; flipud(x(end-(Np-N0)+1:end))];

h = wavelet_filter(cfg.wavelet);

% ---- analysis ----
a = xp;
D = cell(J, 1);
for j = 1:J
    [a, D{j}] = dwt_periodic(a, h);
end

% ---- thresholding ----
sigma_global = median(abs(D{1})) / 0.6745;
Dt  = D;
thr = zeros(J, 1);
for j = 1:J
    switch cfg.wav_rule
        case 'bayes'
            sx = sqrt(max(mean(D{j}.^2) - sigma_global^2, 0));
            if sx > 0
                thr(j) = sigma_global^2 / sx;
            else
                thr(j) = max(abs(D{j}));       % band is pure noise -> remove it
            end
        case 'level'
            sigma  = median(abs(D{j})) / 0.6745;
            thr(j) = sigma * sqrt(2 * log(length(xp)));
        otherwise   % 'global'
            thr(j) = sigma_global * sqrt(2 * log(length(xp)));
    end
    if strcmp(cfg.wav_threshold, 'soft')
        Dt{j} = sign(D{j}) .* max(abs(D{j}) - thr(j), 0);
    else
        Dt{j} = D{j} .* (abs(D{j}) > thr(j));
    end
end

% ---- synthesis ----
y = a;
for j = J:-1:1
    y = idwt_periodic(y, Dt{j}, h);
end
y = y(1:N0);

info.approx   = a;
info.details  = D;
info.details_thr = Dt;
info.thresholds  = thr;
end
