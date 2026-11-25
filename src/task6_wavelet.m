function task6_wavelet(clean, noisy, fs, cfg)
% TASK 6 - Wavelet (DWT) denoising: decomposition, thresholding, variants.
%
%   Uses wavelet_denoise.m (own DWT, implemented from scratch). This task shows
%   the decomposition of the white-noise case and compares threshold rules
%   (BayesShrink / VisuShrink, soft / hard), the number of levels and the
%   wavelet family (Haar, Daubechies, Symlets, Coiflets).

fprintf('\n=== Task 6: Wavelet denoising ===\n');
iw = find(strcmp({noisy.name}, 'white'));
x  = noisy(iw).x;
[y, info] = wavelet_denoise(x, cfg);

%% Fig 14 - detail coefficients before/after thresholding
J = cfg.wav_levels;
fig = new_figure(cfg, 'DWT decomposition', 1100, 900);
for j = 1:J
    subplot(J+1, 1, j);
    d = info.details{j};
    plot(d, 'Color', [0.7 0.7 0.7]); hold on;
    plot(info.details_thr{j}, 'b');
    plot([1 length(d)], info.thresholds(j)*[1 1], 'r--', [1 length(d)], -info.thresholds(j)*[1 1], 'r--');
    hold off; xlim([1 length(d)]); grid on;
    band = fs ./ 2.^[j+1 j];
    ylabel(sprintf('d_%d', j));
    title(sprintf('Detail level %d (%.0f - %.0f Hz): grey = noisy, blue = after %s threshold', ...
          j, band(1), band(2), cfg.wav_threshold));
end
subplot(J+1, 1, J+1);
plot(info.approx, 'k'); xlim([1 length(info.approx)]); grid on;
ylabel(sprintf('a_%d', J));
title(sprintf('Approximation a_%d (0 - %.0f Hz), not thresholded', J, fs/2^(J+1)));
save_figure(fig, cfg, '14_wavelet_decomposition');

%% Fig 15 - time-domain zoom
t = (0:length(x)-1)/fs;
% centre a 300 ms window on the loudest part of the speech
e = filter(ones(round(0.3*fs), 1), 1, clean.^2);
[~, iE] = max(e);
t1 = max(0, t(iE) - 0.3);
seg = t >= t1 & t <= t1 + 0.3;
fig = new_figure(cfg, 'Wavelet time domain', 1100, 450);
plot(t(seg), x(seg), 'Color', [0.8 0.8 0.8]); hold on;
plot(t(seg), clean(seg), 'b', 'LineWidth', 1.2);
plot(t(seg), y(seg), 'r', 'LineWidth', 1); hold off; grid on;
legend('Noisy (white)', 'Clean', 'Wavelet denoised');
xlabel('Time (s)'); ylabel('Amplitude'); title('Wavelet denoising, 300 ms zoom');
save_figure(fig, cfg, '15_wavelet_time_zoom');

%% Variant comparison table
rules = {'bayes', 'global', 'level'};
thr   = {'soft', 'hard'};
levels = 3:7;
fid = fopen(fullfile(cfg.out_dir, 'task6_wavelet_variants.txt'), 'w');
l = sprintf('%-8s %-8s %-6s %6s | %s\n', 'Noise', 'Rule', 'Thr', 'Level', 'dSNR(dB)  LSD(dB)');
fprintf('%s', l); fprintf(fid, '%s', l);
best = struct('dsnr', -Inf);
for k = 1:numel(noisy)
    for r = 1:numel(rules)
        for s = 1:2
            for J = levels
                c = cfg; c.wav_rule = rules{r}; c.wav_threshold = thr{s}; c.wav_levels = J;
                yy = wavelet_denoise(noisy(k).x, c);
                m  = compute_metrics(clean, noisy(k).x, yy, fs, cfg);
                l = sprintf('%-8s %-8s %-6s %6d | %8.2f %8.2f\n', noisy(k).name, ...
                            rules{r}, thr{s}, J, m.dsnr, m.lsd);
                fprintf(fid, '%s', l);
                if k == iw && m.dsnr > best.dsnr
                    best = struct('dsnr', m.dsnr, 'rule', rules{r}, 'thr', thr{s}, 'J', J);
                end
            end
        end
    end
end
l = sprintf('\nBest for white noise: rule=%s, threshold=%s, levels=%d (dSNR = %.2f dB)\n', ...
            best.rule, best.thr, best.J, best.dsnr);
fprintf('%s', l); fprintf(fid, '%s', l);

%% Wavelet family comparison (Haar, Daubechies, Symlets, Coiflets)
fams = cfg.wav_families;
DS = zeros(numel(fams), numel(noisy)); LS = DS;
for w = 1:numel(fams)
    c = cfg; c.wavelet = fams{w};
    for k = 1:numel(noisy)
        yy = wavelet_denoise(noisy(k).x, c);
        m  = compute_metrics(clean, noisy(k).x, yy, fs, cfg);
        DS(w, k) = m.dsnr; LS(w, k) = m.lsd;
    end
end
l = sprintf('\nWavelet family comparison (%s rule, %s threshold, %d levels): dSNR (dB)\n%-8s', ...
            cfg.wav_rule, cfg.wav_threshold, cfg.wav_levels, 'Wavelet');
l = [l sprintf(' %8s', noisy.name) sprintf('\n')];
for w = 1:numel(fams)
    l = [l sprintf('%-8s', fams{w}) sprintf(' %8.2f', DS(w, :)) sprintf('   (%d taps)\n', numel(wavelet_filter(fams{w})))]; %#ok<AGROW>
end
fprintf('%s', l); fprintf(fid, '%s', l);
fclose(fid);

fig = new_figure(cfg, 'Wavelet families', 1000, 450);
bar(DS); grid on;
set(gca, 'XTickLabel', fams);
ylabel('SNR improvement (dB)'); legend({noisy.name}, 'Location', 'eastoutside');
title('Wavelet family comparison');
save_figure(fig, cfg, '15b_wavelet_families');

%% Optional cross-check with the Wavelet Toolbox (wavedec, as in the lab)
try
    if exist('wavedec', 'file') && exist('dwtmode', 'file') && ~exist('OCTAVE_VERSION', 'builtin')
        old = dwtmode('status', 'nodisp');
        dwtmode('per', 'nodisp');
        [~, info2] = wavelet_denoise(x(1:2^14), cfg);
        [cc, ll] = wavedec(x(1:2^14), cfg.wav_levels, cfg.wavelet);
        e_own = sum(cellfun(@(d) sum(d.^2), info2.details)) + sum(info2.approx.^2);
        fprintf('   energy check: own DWT %.6f, wavedec %.6f, signal %.6f\n', ...
                e_own, sum(cc.^2), sum(x(1:2^14).^2));
        dwtmode(old, 'nodisp');
    end
catch
end
fprintf('   full variant table -> results/task6_wavelet_variants.txt\n');
end
