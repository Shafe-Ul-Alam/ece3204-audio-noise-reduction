function plot_filter_set(flts, fs, cfg, ttl_fmt, ttl_arg, fname)
% PLOT_FILTER_SET  Overlay magnitude (full + pass-band zoom) and phase of filters.
fig = new_figure(cfg, 'Filter responses', 1100, 750);
cols = lines(numel(flts));
for k = 1:numel(flts)
    [H, f] = filter_response(flts(k), fs);
    M = 20*log10(abs(H) + 1e-12);

    subplot(2, 2, [1 2]);
    plot(f, M, 'Color', cols(k, :), 'LineWidth', 1.2); hold on;

    subplot(2, 2, 3);
    plot(f, M, 'Color', cols(k, :), 'LineWidth', 1.2); hold on;

    subplot(2, 2, 4);
    plot(f, unwrap(angle(H)), 'Color', cols(k, :), 'LineWidth', 1.2); hold on;
end
names = {flts.name};

subplot(2, 2, [1 2]); hold off; grid on;
ylim([-120 5]); xlim([0 fs/2]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title(sprintf(ttl_fmt, ttl_arg));
legend(names, 'Location', 'southwest');

subplot(2, 2, 3); hold off; grid on;
xlim([cfg.band(1)*0.5, cfg.band(2)*1.2]); ylim([-6 2]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Pass-band detail');

subplot(2, 2, 4); hold off; grid on;
xlim([0 fs/2]);
xlabel('Frequency (Hz)'); ylabel('Phase (rad)');
title('Unwrapped phase');

save_figure(fig, cfg, fname);
end
