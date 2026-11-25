function mr = task4_multirate(clean, fs, cfg)
% TASK 4 - Downsampling / upsampling with anti-aliasing and interpolation filters.
%
% Decimation by M  : y[n] = x[nM]. Any content above fs/(2M) folds back
%                    (aliasing) unless a low-pass anti-aliasing filter with
%                    cut-off pi/M is applied first.
% Interpolation by L: insert L-1 zeros between samples. This creates L-1
%                    spectral images; a low-pass interpolation filter with
%                    cut-off pi/L and gain L removes them.
%
% Part A uses a linear chirp (0 -> fs/2) because aliasing/imaging are then
% obvious in a spectrogram. Part B measures the effect on real speech for
% every factor in cfg.mr_factors.

fprintf('\n=== Task 4: Multirate DSP ===\n');

%% -------- Part A: chirp demonstration (M = L = 4) --------
M  = 4;
tc = (0:2*fs-1)' / fs;
chirp_x = 0.5 * cos(2*pi*(0.5*(fs/2)/2) * tc.^2);      % 0 -> fs/2 in 2 s

[aa, ip] = mr_filters(M, fs, cfg);
d_naive = chirp_x(1:M:end);                            % no anti-aliasing
d_aa    = apply_filter(aa, chirp_x); d_aa = d_aa(1:M:end);
u_zero  = zero_insert(d_aa, M);                         % imaging
u_interp = apply_filter(ip, u_zero);                    % images removed

fig = new_figure(cfg, 'Multirate chirp demo', 1150, 750);
subplot(2, 2, 1); plot_spectrogram(d_naive, fs/M, sprintf('Downsample x%d WITHOUT anti-alias filter (aliasing)', M));
subplot(2, 2, 2); plot_spectrogram(d_aa,    fs/M, sprintf('Downsample x%d WITH anti-alias filter', M));
subplot(2, 2, 3); plot_spectrogram(u_zero,  fs,   sprintf('Upsample x%d zero insertion only (imaging)', M));
subplot(2, 2, 4); plot_spectrogram(u_interp, fs,  sprintf('Upsample x%d WITH interpolation filter', M));
save_figure(fig, cfg, '08_multirate_chirp_aliasing_imaging');

%% -------- Part B: speech, all factors --------
F = cfg.mr_factors;
mr.factor = F(:);
mr.snr_naive = zeros(numel(F), 1);  % no anti-alias filter, proper interpolation
mr.snr_aa    = zeros(numel(F), 1);  % anti-alias + interpolation (correct system)
mr.snr_img   = zeros(numel(F), 1);  % anti-alias, zero insertion only
mr.lsd_naive = zeros(numel(F), 1);
mr.lsd_aa    = zeros(numel(F), 1);
mr.bw_loss   = zeros(numel(F), 1);  % % of speech energy above new Nyquist

[Pc, fw] = psd_welch(clean, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);

for i = 1:numel(F)
    M = F(i);
    [aa, ip] = mr_filters(M, fs, cfg);

    dn = clean(1:M:end);                         % naive decimation
    da = apply_filter(aa, clean); da = da(1:M:end);

    y_naive = apply_filter(ip, zero_insert(dn, M));
    y_aa    = apply_filter(ip, zero_insert(da, M));
    y_img   = M * zero_insert(da, M);            % gain M, but no image removal

    n = min(length(clean), length(y_aa));
    m1 = compute_metrics(clean(1:n), clean(1:n), y_naive(1:n), fs, cfg);
    m2 = compute_metrics(clean(1:n), clean(1:n), y_aa(1:n),    fs, cfg);
    m3 = compute_metrics(clean(1:n), clean(1:n), y_img(1:n),   fs, cfg);
    mr.snr_naive(i) = m1.snr_out;  mr.lsd_naive(i) = m1.lsd;
    mr.snr_aa(i)    = m2.snr_out;  mr.lsd_aa(i)    = m2.lsd;
    mr.snr_img(i)   = m3.snr_out;
    mr.bw_loss(i)   = 100 * sum(Pc(fw > fs/(2*M))) / sum(Pc);

    write_wav(da,      fs/M, cfg, sprintf('multirate_down%d_antialias', M));
    write_wav(dn,      fs/M, cfg, sprintf('multirate_down%d_naive', M));
    write_wav(y_aa,    fs,   cfg, sprintf('multirate_roundtrip%d_proper', M));
    write_wav(y_naive, fs,   cfg, sprintf('multirate_roundtrip%d_aliased', M));
    write_wav(y_img,   fs,   cfg, sprintf('multirate_roundtrip%d_imaging', M));

    if M == max(F)          % spectra for the largest factor
        fig = new_figure(cfg, 'Multirate speech spectra', 1150, 750);
        [X0, f0] = fft_db(clean, fs);
        [Xn, fn] = fft_db(dn, fs/M);
        [Xa, fa] = fft_db(da, fs/M);
        [Xi, fi] = fft_db(y_img, fs);
        [Xu, fu] = fft_db(y_aa, fs);
        subplot(2, 2, 1); plot(f0, X0); hold on; plot(fs/(2*M)*[1 1], [-150 0], 'r--'); hold off;
        grid on; xlim([0 fs/2]); title('Original speech (red = new Nyquist)'); xlabel('Hz'); ylabel('dB');
        subplot(2, 2, 2); plot(fn, Xn, 'r'); hold on; plot(fa, Xa, 'b'); hold off; grid on;
        xlim([0 fs/(2*M)]); legend('no anti-alias', 'with anti-alias');
        title(sprintf('Downsampled x%d (fs = %d Hz)', M, fs/M)); xlabel('Hz'); ylabel('dB');
        subplot(2, 2, 3); plot(fi, Xi, 'm'); grid on; xlim([0 fs/2]);
        title(sprintf('Upsampled x%d, zero insertion: spectral images', M)); xlabel('Hz'); ylabel('dB');
        subplot(2, 2, 4); plot(fu, Xu, 'b'); grid on; xlim([0 fs/2]);
        title(sprintf('Upsampled x%d with interpolation filter', M)); xlabel('Hz'); ylabel('dB');
        save_figure(fig, cfg, '09_multirate_speech_spectra');
    end
end

% bar chart of round-trip quality
fig = new_figure(cfg, 'Multirate SNR', 800, 450);
bar([mr.snr_aa mr.snr_naive mr.snr_img]); grid on;
set(gca, 'XTickLabel', arrayfun(@(m) sprintf('M = L = %d', m), F, 'UniformOutput', false));
ylabel('Round-trip SNR vs original (dB)');
legend('anti-alias + interpolation', 'no anti-alias filter', 'no interpolation filter', 'Location', 'northeast');
title('Down/up-sampling round trip on clean speech');
save_figure(fig, cfg, '10_multirate_snr');

% text table
fid = fopen(fullfile(cfg.out_dir, 'task4_multirate.txt'), 'w');
hdr = sprintf('%6s %10s %14s %14s %14s %12s %12s\n', 'Factor', 'New fs', ...
              'SNR proper', 'SNR aliased', 'SNR imaging', 'LSD proper', 'Energy lost');
fprintf('%s', hdr); fprintf(fid, '%s', hdr);
for i = 1:numel(F)
    l = sprintf('%6d %10d %14.2f %14.2f %14.2f %12.2f %11.2f%%\n', F(i), fs/F(i), ...
                mr.snr_aa(i), mr.snr_naive(i), mr.snr_img(i), mr.lsd_aa(i), mr.bw_loss(i));
    fprintf('%s', l); fprintf(fid, '%s', l);
end
fclose(fid);
end

% --------------------------------------------------------------------------
function [aa, ip] = mr_filters(M, fs, cfg)
% Anti-aliasing (cut-off 0.9*pi/M) and interpolation (same cut-off, gain M) FIR filters.
% apply_filter runs the filter forward AND backward (zero phase), so the gain
% is applied twice: sqrt(M) per pass gives the required overall gain of M.
b = fir_window_design(cfg.mr_fir_order + 1, 0.9 * (fs/2) / M, fs, 'lowpass', 'hamming');
aa = struct('name', 'antialias', 'type', 'fir', 'b', b, 'sos', [], 'g', 1);
ip = struct('name', 'interp',    'type', 'fir', 'b', sqrt(M) * b, 'sos', [], 'g', 1);
end

function y = zero_insert(x, L)
y = zeros(length(x) * L, 1);
y(1:L:end) = x;
end
