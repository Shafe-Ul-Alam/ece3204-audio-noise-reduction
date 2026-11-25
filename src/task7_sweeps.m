function task7_sweeps(clean, noisy, fs, cfg)
% TASK 7 (part 2) - Effect of filter order and cut-off (design challenge 2).
%
% On the mixed-noise signal:
%   (a) FIR-Hamming order sweep      (cfg.sweep_fir_orders)
%   (b) IIR-Butterworth order sweep  (cfg.sweep_iir_orders)
%   (c) Upper cut-off sweep for both (cfg.sweep_upper_fc)
% Reports SNR improvement, log-spectral distance and runtime.

fprintf('\n=== Task 7b: Order and cut-off sweeps (mixed noise) ===\n');
x = noisy(end).x;
fid = fopen(fullfile(cfg.out_dir, 'task7_sweeps.txt'), 'w');

% (a) FIR order
fo = cfg.sweep_fir_orders;
[dF, lF, tF] = deal(zeros(size(fo)));
for i = 1:numel(fo)
    flt = make_fir(fo(i), cfg.band, fs);
    [dF(i), lF(i), tF(i)] = run_one(flt, clean, x, fs, cfg);
end
% (b) IIR order
io = cfg.sweep_iir_orders;
[dI, lI, tI] = deal(zeros(size(io)));
for i = 1:numel(io)
    flt = make_butter(io(i), cfg.band, fs);
    [dI(i), lI(i), tI(i)] = run_one(flt, clean, x, fs, cfg);
end
% (c) cut-off
fc = cfg.sweep_upper_fc;
[dC1, lC1, dC2, lC2] = deal(zeros(size(fc)));
for i = 1:numel(fc)
    b = [cfg.band(1) fc(i)];
    [dC1(i), lC1(i)] = run_one(make_fir(cfg.fir_order, b, fs), clean, x, fs, cfg);
    [dC2(i), lC2(i)] = run_one(make_butter(cfg.iir_order, b, fs), clean, x, fs, cfg);
end

% ---- text output ----
p = @(varargin) fprintf(fid, varargin{:});
p('(a) FIR-Hamming order sweep\n%8s %10s %10s %10s\n', 'Order', 'dSNR', 'LSD', 'ms');
p('%8d %10.2f %10.2f %10.2f\n', [fo; dF; lF; tF]);
p('\n(b) IIR-Butterworth order sweep (prototype order; band-pass order = 2x)\n%8s %10s %10s %10s\n', 'Order', 'dSNR', 'LSD', 'ms');
p('%8d %10.2f %10.2f %10.2f\n', [io; dI; lI; tI]);
p('\n(c) Upper cut-off sweep (lower edge %d Hz)\n%8s %12s %12s %12s %12s\n', cfg.band(1), 'fc(Hz)', 'FIR dSNR', 'FIR LSD', 'IIR dSNR', 'IIR LSD');
p('%8d %12.2f %12.2f %12.2f %12.2f\n', [fc; dC1; lC1; dC2; lC2]);
fclose(fid);
type(fullfile(cfg.out_dir, 'task7_sweeps.txt'));

% ---- Fig 20 ----
fig = new_figure(cfg, 'Sweeps', 1200, 750);
subplot(2, 2, 1);
plot(fo, dF, 'o-', 'LineWidth', 1.4); grid on;
xlabel('FIR order'); ylabel('dSNR (dB)'); title('FIR-Hamming: order vs noise reduction');
subplot(2, 2, 2);
plot(io, dI, 's-', 'LineWidth', 1.4); grid on;
xlabel('Butterworth prototype order'); ylabel('dSNR (dB)'); title('Butterworth: order vs noise reduction');
subplot(2, 2, 3);
plot(tF, dF, 'o-', 'LineWidth', 1.4); hold on; plot(tI, dI, 's-', 'LineWidth', 1.4);
for i = 1:numel(fo), text(tF(i), dF(i) - 0.25, sprintf('N=%d', fo(i)), 'FontSize', 8); end
for i = 1:numel(io), text(tI(i), dI(i) + 0.25, sprintf('n=%d', io(i)), 'FontSize', 8); end
hold off; grid on;
xlabel('Runtime (ms)'); ylabel('dSNR (dB)'); legend('FIR-Hamming', 'Butterworth', 'Location', 'southeast');
title('Noise reduction vs computational cost');
subplot(2, 2, 4);
plot(fc, dC1, 'o-', fc, dC2, 's-', 'LineWidth', 1.4); hold on;
plot(fc, lC1, 'o--', fc, lC2, 's--', 'LineWidth', 1.2); hold off; grid on;
xlabel('Upper cut-off (Hz)'); ylabel('dB');
legend('FIR dSNR', 'IIR dSNR', 'FIR LSD', 'IIR LSD', 'Location', 'east');
title('Cut-off: noise reduction vs distortion');
save_figure(fig, cfg, '20_order_cutoff_sweeps');
end

% --------------------------------------------------------------------------
function [d, l, t] = run_one(flt, clean, x, fs, cfg)
tt = zeros(cfg.timing_runs, 1);
for r = 1:cfg.timing_runs
    t0 = tic; y = apply_filter(flt, x); tt(r) = toc(t0);
end
m = compute_metrics(clean, x, y, fs, cfg);
d = m.dsnr; l = m.lsd; t = 1000*mean(tt);
end

function flt = make_fir(N, band, fs)
flt = struct('name', 'fir', 'type', 'fir', 'b', fir_window_design(N + 1, band, fs, 'bandpass', 'hamming'), 'sos', [], 'g', 1);
end

function flt = make_butter(n, band, fs)
[z, p, k] = butter(n, band/(fs/2));
[sos, g] = zp2sos(z, p, k);
flt = struct('name', 'butter', 'type', 'iir', 'b', [], 'sos', sos, 'g', g);
end
