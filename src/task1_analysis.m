function task1_analysis(clean, noisy, fs, cfg)
% TASK 1 - Load WAV, plot waveform, inspect FFT spectrum (DFT and FFT).
%
% Also verifies a DFT computed from its equation and an own radix-2 DIT FFT
% against fft(), compares their cost, and shows FFT-based (fast) convolution.
%
% Outputs
%   fig 01 : time-domain waveforms (clean + every noisy version)
%   fig 02 : FFT magnitude spectra, noisy vs clean
%   fig 03 : spectrograms of clean and mixed-noise signals
%   results/task1_analysis.txt : speech band + detected noise components

fprintf('\n=== Task 1: Time-domain and FFT analysis ===\n');
t = (0:length(clean)-1) / fs;
nN = numel(noisy);

%% Fig 01 - waveforms
fig = new_figure(cfg, 'Waveforms', 1000, 800);
subplot(nN+1, 1, 1);
plot(t, clean, 'b'); grid on; ylim([-1 1]);
title('Clean reference speech'); ylabel('Amplitude');
for k = 1:nN
    subplot(nN+1, 1, k+1);
    plot(t, noisy(k).x, 'r'); grid on; ylim([-1 1]);
    title(sprintf('Noisy: %s noise (input SNR = %d dB)', noisy(k).name, noisy(k).snr_target));
    ylabel('Amplitude');
end
xlabel('Time (s)');
save_figure(fig, cfg, '01_waveforms');

%% Fig 02 - FFT spectra
[Xc, f] = fft_db(clean, fs);
fig = new_figure(cfg, 'FFT spectra', 1100, 750);
for k = 1:nN
    subplot(2, 2, k);
    Xn = fft_db(noisy(k).x, fs);
    plot(f/1000, Xn, 'r'); hold on;
    plot(f/1000, Xc, 'b'); hold off; grid on;
    xlim([0 fs/2000]); ylim([max(Xc)-110, max([Xc; Xn])+5]);
    xlabel('Frequency (kHz)'); ylabel('|X(f)| (dB)');
    title(sprintf('FFT: %s noise', noisy(k).name));
    legend('Noisy', 'Clean', 'Location', 'northeast');
end
save_figure(fig, cfg, '02_fft_spectra');

%% Fig 03 - spectrograms
fig = new_figure(cfg, 'Spectrograms', 1100, 450);
subplot(1, 2, 1); plot_spectrogram(clean, fs, 'Clean speech');
subplot(1, 2, 2); plot_spectrogram(noisy(end).x, fs, sprintf('Noisy (%s)', noisy(end).name));
save_figure(fig, cfg, '03_spectrograms');

%% Numerical analysis: where is the speech, where is the noise?
fid = fopen(fullfile(cfg.out_dir, 'task1_analysis.txt'), 'w');
out = @(varargin) fprintf_both(fid, varargin{:});

% speech band from cumulative energy of the clean spectrum
[Pc, fw] = psd_welch(clean, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
cE = cumsum(Pc) / sum(Pc);
f_lo = fw(find(cE >= 0.01, 1)); f_hi = fw(find(cE >= 0.99, 1));
out('Signal length: %.2f s, fs = %d Hz, N = %d samples\n', length(clean)/fs, fs, length(clean));
out('98%% of clean speech energy lies between %.0f Hz and %.0f Hz\n', f_lo, f_hi);
out('Chosen pass-band for filters: %d - %d Hz\n\n', cfg.band(1), cfg.band(2));

% noise components: frequencies where the noisy PSD exceeds the clean PSD
for k = 1:nN
    Pn = psd_welch(noisy(k).x, fs, cfg.welch_win, cfg.welch_overlap, cfg.nfft);
    excess = 10*log10(Pn + eps) - 10*log10(Pc + eps);
    % local maxima of the excess that are > 15 dB
    pk = find(excess(2:end-1) > excess(1:end-2) & excess(2:end-1) >= excess(3:end) ...
              & excess(2:end-1) > 15) + 1;
    [~, order] = sort(excess(pk), 'descend');
    pk = pk(order(1:min(6, numel(order))));
    frac_low  = sum(Pn(fw < cfg.band(1)))  / sum(Pn);
    frac_high = sum(Pn(fw > cfg.band(2)))  / sum(Pn);
    out('[%s] strongest noise peaks (Hz): %s\n', noisy(k).name, mat2str(round(sort(fw(pk)'))));
    out('        noisy power below %d Hz: %.1f %%, above %d Hz: %.1f %%\n', ...
        cfg.band(1), 100*frac_low, cfg.band(2), 100*frac_high);
end

%% DFT vs FFT (lab experiment: DFT by equation, DIT FFT, circular convolution)
out('\n--- DFT / FFT verification on a 1024-sample speech frame ---\n');
[~, i0] = max(abs(clean));
fr = clean(max(1, i0-512) + (0:1023));
X1 = dft_direct(fr); X2 = fft_dit(fr); X3 = fft(fr);
out('max |DFT(equation) - fft()|      = %.2e\n', max(abs(X1 - X3)));
out('max |FFT(own radix-2 DIT) - fft()| = %.2e\n', max(abs(X2 - X3)));
out('Parseval: sum|x|^2 = %.6f, (1/N) sum|X|^2 = %.6f\n', sum(fr.^2), sum(abs(X3).^2)/1024);

Ns = 2.^(6:11);
[tD, tF, tB] = deal(zeros(size(Ns)));
for i = 1:numel(Ns)
    v = randn(Ns(i), 1);
    t0 = tic; for r = 1:3, dft_direct(v); end; tD(i) = toc(t0)/3;
    t0 = tic; for r = 1:3, fft_dit(v);    end; tF(i) = toc(t0)/3;
    t0 = tic; for r = 1:50, fft(v);       end; tB(i) = toc(t0)/50;
end
out('\n%8s %14s %14s %14s %12s %12s\n', 'N', 'DFT eq. (ms)', 'own FFT (ms)', 'fft() (ms)', 'N^2 mult', '(N/2)log2N');
for i = 1:numel(Ns)
    out('%8d %14.3f %14.3f %14.4f %12d %12d\n', Ns(i), 1e3*tD(i), 1e3*tF(i), 1e3*tB(i), ...
        Ns(i)^2, Ns(i)/2*log2(Ns(i)));
end

fig = new_figure(cfg, 'DFT vs FFT', 1100, 450);
subplot(1, 2, 1);
loglog(Ns, 1e3*tD, 'ro-', Ns, 1e3*tF, 'bs-', Ns, 1e3*tB, 'k^-', 'LineWidth', 1.3); grid on;
xlabel('N'); ylabel('Time (ms)'); legend('DFT by equation', 'Own radix-2 DIT FFT', 'Built-in fft()', 'Location', 'northwest');
title('Measured computation time');
subplot(1, 2, 2);
loglog(Ns, Ns.^2, 'ro-', Ns, Ns/2.*log2(Ns), 'bs-', 'LineWidth', 1.3); grid on;
xlabel('N'); ylabel('Complex multiplications'); legend('DFT: N^2', 'FFT: (N/2) log_2 N', 'Location', 'northwest');
title('Theoretical cost');
save_figure(fig, cfg, '03b_dft_vs_fft');

% FFT-based (fast) convolution = FIR filtering; circular wrap-around if not zero-padded
h  = fir_window_design(cfg.fir_order + 1, cfg.band, fs, 'bandpass', 'hamming');
xs = noisy(end).x(1:8192);
y_lin  = filter(h, 1, xs);                                   % direct linear convolution
L      = 2^nextpow2(length(xs) + length(h) - 1);
y_fft  = real(ifft(fft(xs, L) .* fft(h(:), L)));             % zero-padded -> linear
y_circ = real(ifft(fft(xs) .* fft(h(:), length(xs))));       % no padding -> circular
out('\nFast convolution via FFT (8192 samples, %d-tap FIR):\n', length(h));
out('  zero-padded FFT vs filter():   max error = %.2e  (linear convolution)\n', ...
    max(abs(y_fft(1:length(xs)) - y_lin)));
out('  no zero-padding vs filter():   max error = %.2e  in first %d samples (circular wrap-around)\n', ...
    max(abs(y_circ(1:length(h)) - y_lin(1:length(h)))), length(h));
fclose(fid);
end

function fprintf_both(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
