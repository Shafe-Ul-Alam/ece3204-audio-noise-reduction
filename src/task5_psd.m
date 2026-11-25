function task5_psd(clean, noisy, fs, cfg)
% TASK 5 - Power spectral density: Periodogram vs Welch.
%
%   Periodogram : single |FFT|^2 of the whole record. Resolution fs/N is very
%                 fine, but the variance does not decrease with N (each bin is
%                 a chi-square variable with 2 degrees of freedom).
%   Welch       : average of windowed periodograms of overlapping segments.
%                 Variance drops roughly as 1/K (K = number of segments) but
%                 resolution becomes ~fs/segment_length.

fprintf('\n=== Task 5: PSD estimation (Periodogram vs Welch) ===\n');
xm = noisy(end).x;                                    % mixed-noise signal

%% Fig 11 - periodogram vs Welch of the mixed-noise signal
[Pp, fp] = psd_periodogram(xm, fs);
[Pw, fw] = psd_welch(xm, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
fig = new_figure(cfg, 'Periodogram vs Welch', 1100, 500);
plot(fp, 10*log10(Pp + eps), 'Color', [0.7 0.7 0.7]); hold on;
plot(fw, 10*log10(Pw + eps), 'b', 'LineWidth', 1.6); hold off; grid on;
xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
legend('Periodogram (rectangular, whole record)', ...
       sprintf('Welch (Hamming %d, %d%% overlap)', cfg.welch_win, round(100*cfg.welch_overlap)));
title(sprintf('PSD of speech + %s noise', noisy(end).name));
save_figure(fig, cfg, '11_psd_periodogram_vs_welch');

%% Fig 12 - variance / resolution trade-off on pure white noise
rng(cfg.seed + 1);
wn = randn(length(clean), 1);           % true PSD = 2/fs (one-sided, unit variance)
true_db = 10*log10(2/fs);
segs = [128 512 2048];
fig = new_figure(cfg, 'Welch trade-off', 1100, 500);
[Pp, fp] = psd_periodogram(wn, fs);
plot(fp, 10*log10(Pp), 'Color', [0.75 0.75 0.75]); hold on;
cols = lines(numel(segs));
fid = fopen(fullfile(cfg.out_dir, 'task5_psd.txt'), 'w');
l = sprintf('White-noise test (true PSD = %.2f dB/Hz)\n%-28s %10s %14s %14s\n', ...
            true_db, 'Estimator', 'Segments', 'Resolution(Hz)', 'Std dev (dB)');
fprintf('%s', l); fprintf(fid, '%s', l);
inner = fp > 100 & fp < fs/2 - 100;
l = sprintf('%-28s %10d %14.2f %14.2f\n', 'Periodogram', 1, fs/length(wn), std(10*log10(Pp(inner))));
fprintf('%s', l); fprintf(fid, '%s', l);
leg = {'Periodogram'};
for i = 1:numel(segs)
    [P, f] = psd_welch(wn, fs, segs(i), 0.5, segs(i));
    plot(f, 10*log10(P), 'Color', cols(i, :), 'LineWidth', 1.4);
    K = floor((length(wn) - segs(i)) / (segs(i)/2)) + 1;
    in = f > 100 & f < fs/2 - 100;
    l = sprintf('%-28s %10d %14.2f %14.2f\n', sprintf('Welch, segment %d', segs(i)), ...
                K, 1.36*fs/segs(i), std(10*log10(P(in))));
    fprintf('%s', l); fprintf(fid, '%s', l);
    leg{end+1} = sprintf('Welch L = %d (K = %d)', segs(i), K); %#ok<AGROW>
end
plot([0 fs/2], true_db*[1 1], 'k--', 'LineWidth', 1.5);
hold off; grid on; xlim([0 fs/2]); ylim([true_db-30 true_db+12]);
leg{end+1} = 'True PSD';
legend(leg, 'Location', 'southeast');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Variance vs resolution: white Gaussian noise');
save_figure(fig, cfg, '12_psd_variance_tradeoff');
fprintf(fid, '\nResolution of Welch with a Hamming window ~ 1.36 fs / L (main-lobe 3 dB width).\n');
fclose(fid);

%% Fig 12b - parametric AR PSD (Yule-Walker, Burg) vs Welch (lab experiment)
fid = fopen(fullfile(cfg.out_dir, 'task5_psd.txt'), 'a');
fig = new_figure(cfg, 'AR PSD', 1150, 500);
[Pw, fw] = psd_welch(xm, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
plot(fw, 10*log10(Pw), 'Color', [0.6 0.6 0.6], 'LineWidth', 2); hold on;
leg = {'Welch (non-parametric)'};
sty = {'b', 'r', 'b--', 'r--'};
i = 0;
fprintf(fid, '\nAR (parametric) PSD of the mixed-noise signal\n%-14s %6s %18s\n', 'Method', 'Order', 'Pred. error var');
for p = cfg.ar_orders
    for meth = {'yule', 'burg'}
        i = i + 1;
        [Pa, fa, ~, E] = ar_psd_estimate(xm, p, cfg.nfft, fs, meth{1});
        plot(fa, 10*log10(Pa), sty{i}, 'LineWidth', 1.2);
        leg{end+1} = sprintf('%s, order %d', strrep(strrep(meth{1}, 'yule', 'Yule-Walker'), 'burg', 'Burg'), p); %#ok<AGROW>
        fprintf(fid, '%-14s %6d %18.3e\n', meth{1}, p, E);
    end
end
hold off; grid on; xlim([0 fs/2]);
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)'); legend(leg, 'Location', 'southwest');
title('Parametric AR PSD vs Welch (speech + mixed noise)');
save_figure(fig, cfg, '12b_psd_ar_yulewalker_burg');

% optional cross-check against the toolbox functions used in the lab
try
    if ~exist('OCTAVE_VERSION', 'builtin')
        [p1, ~] = periodogram(xm, [], length(xm), fs);
        [p2, ~] = pwelch(xm, hamming(cfg.welch_win), round(cfg.welch_overlap*cfg.welch_win), cfg.nfft, fs);
        [p3, ~] = pburg(xm, cfg.ar_orders(end), cfg.nfft, fs);
        [q1, ~] = psd_periodogram(xm, fs, length(xm));
        [q3, ~] = ar_psd_estimate(xm, cfg.ar_orders(end), cfg.nfft, fs, 'burg');
        d = @(a, b) median(abs(10*log10(a(:) + eps) - 10*log10(b(:) + eps)));
        fprintf(fid, ['\nCheck vs toolbox (median dB difference): periodogram %.3f, ' ...
                      'pwelch %.3f, pburg %.3f\n'], d(p1, q1), d(p2, Pw), d(p3, q3));
    end
catch err
    fprintf(fid, '\nToolbox cross-check skipped: %s\n', err.message);
end
fclose(fid);

%% Fig 13 - Welch PSD of every noisy signal vs clean
[Pc, fc] = psd_welch(clean, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
fig = new_figure(cfg, 'Welch PSDs', 1100, 750);
for k = 1:numel(noisy)
    subplot(2, 2, k);
    Pn = psd_welch(noisy(k).x, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
    plot(fc, 10*log10(Pn + eps), 'r'); hold on;
    plot(fc, 10*log10(Pc + eps), 'b'); hold off; grid on;
    xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title(sprintf('Welch PSD: %s noise', noisy(k).name));
    legend('Noisy', 'Clean', 'Location', 'northeast');
end
save_figure(fig, cfg, '13_psd_welch_all_noises');
end
