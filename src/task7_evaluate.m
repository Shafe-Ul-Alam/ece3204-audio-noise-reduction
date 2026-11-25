function R = task7_evaluate(clean, noisy, fs, filters, cfg)
% TASK 7 - Apply every method to every noisy signal and compare them.
%
% Methods : baseline (no processing), 3 FIR + 3 IIR band-pass filters, wavelet denoising, a hybrid
%           (FIR-Hamming band-pass followed by wavelet) and the adaptive NLMS
%           hum canceller.
% Metrics : input/output SNR, SNR improvement, segmental SNR, log-spectral
%           distance (distortion), correlation, runtime.
% Outputs : results/results_table.csv, results/summary.txt, figures 16-19,
%           processed WAV files in results/audio.

fprintf('\n=== Task 7: Evaluation and comparison ===\n');

% ---- method list (name + function handle) ----
methods = {'No-processing', @(x) x};       % baseline: the noisy input itself
for k = 1:numel(filters)
    methods(end+1, :) = {filters(k).name, @(x) apply_filter(filters(k), x)}; %#ok<AGROW>
end
iHam = find(strcmp({filters.name}, 'FIR-Hamming'));
methods(end+1, :) = {['Wavelet-' cfg.wavelet], @(x) wavelet_denoise(x, cfg)};
methods(end+1, :) = {'Hybrid-FIR+Wavelet', @(x) wavelet_denoise(apply_filter(filters(iHam), x), cfg)};
methods(end+1, :) = {'Adaptive-NLMS',    @(x) nlms_hum_canceller(x, fs, cfg)};
nM = size(methods, 1); nN = numel(noisy);

% ---- run ----
R = struct('noise', {}, 'method', {}, 'snr_in', {}, 'snr_out', {}, 'dsnr', {}, ...
           'segsnr', {}, 'lsd', {}, 'corr', {}, 'time_ms', {});
out_sig = cell(nN, nM);
write_wav(clean, fs, cfg, '00_clean_reference');
for n = 1:nN
    x = noisy(n).x;
    write_wav(x, fs, cfg, sprintf('%s_0_noisy', noisy(n).name));
    for m = 1:nM
        runs = cfg.timing_runs;
        if strcmp(methods{m, 1}, 'Adaptive-NLMS'), runs = 1; end
        tt = zeros(runs, 1);
        for r = 1:runs
            t0 = tic; y = methods{m, 2}(x); tt(r) = toc(t0);
        end
        y = y(:);
        q = compute_metrics(clean, x, y, fs, cfg);
        R(end+1) = struct('noise', noisy(n).name, 'method', methods{m, 1}, ... %#ok<AGROW>
            'snr_in', q.snr_in, 'snr_out', q.snr_out, 'dsnr', q.dsnr, 'segsnr', q.segsnr, ...
            'lsd', q.lsd, 'corr', q.corr, 'time_ms', 1000*mean(tt));
        out_sig{n, m} = y;
        if m > 1
            write_wav(y, fs, cfg, sprintf('%s_%s', noisy(n).name, methods{m, 1}));
        end
        fprintf('   %-6s %-20s dSNR = %6.2f dB  LSD = %5.2f dB  t = %7.1f ms\n', ...
                noisy(n).name, methods{m, 1}, q.dsnr, q.lsd, 1000*mean(tt));
    end
end

% ---- CSV table ----
fid = fopen(fullfile(cfg.out_dir, 'results_table.csv'), 'w');
fprintf(fid, 'noise,method,snr_in_dB,snr_out_dB,snr_improvement_dB,segsnr_dB,lsd_dB,correlation,runtime_ms\n');
for i = 1:numel(R)
    fprintf(fid, '%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.2f\n', R(i).noise, R(i).method, ...
        R(i).snr_in, R(i).snr_out, R(i).dsnr, R(i).segsnr, R(i).lsd, R(i).corr, R(i).time_ms);
end
fclose(fid);
fprintf('   table -> results/results_table.csv\n');

% ---- matrices for plotting ----
D = reshape([R.dsnr],    nM, nN);   % methods x noises
L = reshape([R.lsd],     nM, nN);
T = reshape([R.time_ms], nM, nN);
mnames = methods(:, 1)'; nnames = {noisy.name};

fig = new_figure(cfg, 'SNR improvement', 1200, 500);
bar(D'); grid on;
set(gca, 'XTickLabel', nnames);
ylabel('SNR improvement (dB)'); title('SNR improvement by method and noise type');
legend(mnames, 'Location', 'eastoutside');
save_figure(fig, cfg, '16_snr_improvement');

fig = new_figure(cfg, 'Distortion', 1200, 500);
bar(L'); grid on;
set(gca, 'XTickLabel', nnames);
ylabel('Log-spectral distance (dB)  - lower is better');
title('Spectral distortion of the recovered speech');
legend(mnames, 'Location', 'eastoutside');
save_figure(fig, cfg, '17_lsd_distortion');

fig = new_figure(cfg, 'Runtime', 1000, 450);
barh(mean(T, 2)); grid on;
set(gca, 'YTick', 1:nM, 'YTickLabel', mnames);
xlabel(sprintf('Mean runtime (ms) for %.1f s of audio', length(clean)/fs));
title('Computational cost');
save_figure(fig, cfg, '18_runtime');

% ---- best method per noise (and balanced score) ----
fid = fopen(fullfile(cfg.out_dir, 'summary.txt'), 'w');
pr(fid, '\n================ SUMMARY ================\n');
for n = 1:nN
    c = 2:nM;                                   % exclude the baseline
    [~, iS] = max(D(c, n)); iS = c(iS);
    [~, iL] = min(L(c, n)); iL = c(iL);
    % balanced score: rank in dSNR + rank in LSD + 0.5 * rank in runtime
    score = rankvec(-D(c, n)) + rankvec(L(c, n)) + 0.5 * rankvec(T(c, n));
    [~, iB] = min(score); iB = c(iB);
    pr(fid, '%-6s : noisy LSD = %.2f dB\n', nnames{n}, L(1, n));
    pr(fid, '%-6s : best SNR gain = %-20s (%.2f dB) | least distortion = %-20s (%.2f dB) | best balance = %s\n', ...
       nnames{n}, mnames{iS}, D(iS, n), mnames{iL}, L(iL, n), mnames{iB});
    best_idx(n) = iB; %#ok<AGROW>
end
fclose(fid);

% ---- Fig 19: mixed-noise case, PSD + waveform of the best method ----
n = nN; m = best_idx(n);
[Pc, f] = psd_welch(clean, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
Pn = psd_welch(noisy(n).x, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
Py = psd_welch(out_sig{n, m}, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
t = (0:length(clean)-1)/fs;
fig = new_figure(cfg, 'Best result', 1150, 750);
subplot(2, 1, 1);
plot(f, 10*log10(Pn), 'r', f, 10*log10(Pc), 'b', f, 10*log10(Py), 'k', 'LineWidth', 1.1);
grid on; xlim([0 fs/2]); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
legend('Noisy', 'Clean', mnames{m}); title(sprintf('Welch PSD, %s noise', nnames{n}));
subplot(2, 1, 2);
e = filter(ones(round(0.3*fs), 1), 1, clean.^2);
[~, iE] = max(e); t1 = max(0, t(iE) - 0.3);
seg = t >= t1 & t <= t1 + 0.3;                 % 300 ms around the loudest part
plot(t(seg), noisy(n).x(seg), 'Color', [0.8 0.8 0.8]); hold on;
plot(t(seg), clean(seg), 'b', 'LineWidth', 1.3);
plot(t(seg), out_sig{n, m}(seg), 'r'); hold off; grid on;
legend('Noisy', 'Clean', mnames{m}); xlabel('Time (s)'); ylabel('Amplitude');
title('Waveforms (300 ms zoom)');
save_figure(fig, cfg, '19_best_method_mixed_noise');
end

function r = rankvec(v)
[~, o] = sort(v(:));
r = zeros(numel(v), 1);
r(o) = 1:numel(v);
end

function pr(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
