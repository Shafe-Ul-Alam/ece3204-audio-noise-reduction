function iir = task3_iir_design(fs, cfg)
% TASK 3 - IIR band-pass filters: Butterworth, Chebyshev-I and Elliptic.
%
% All three are designed by the bilinear transform of an analog prototype
% (MATLAB pre-warps the band edges automatically):
%   Butterworth : maximally flat pass-band, slowest roll-off
%   Chebyshev-I : equiripple pass-band (cfg.cheby_Rp dB), faster roll-off
%   Elliptic    : ripple in pass- and stop-band, sharpest roll-off for a
%                 given order (cfg.ellip_Rs dB stop-band)
% The bilinear-transform vs impulse-invariance comparison at the end shows
% why the bilinear transform (with pre-warping) is used.
% A band-pass design of prototype order n has order 2n. Filters are stored
% as second-order sections (SOS) because high-order direct-form IIR
% coefficients are numerically sensitive and can become unstable.

fprintf('\n=== Task 3: IIR filter design (bilinear transform) ===\n');
n  = cfg.iir_order;
Wn = cfg.band / (fs/2);

[z{1}, p{1}, k{1}] = butter(n, Wn);
[z{2}, p{2}, k{2}] = cheby1(n, cfg.cheby_Rp, Wn);
[z{3}, p{3}, k{3}] = ellip(n, cfg.cheby_Rp, cfg.ellip_Rs, Wn);
names = {'IIR-Butterworth', 'IIR-Chebyshev1', 'IIR-Elliptic'};

for i = 1:3
    [sos, g] = zp2sos(z{i}, p{i}, k{i});
    iir(i).name   = names{i};          %#ok<AGROW>
    iir(i).type   = 'iir';             %#ok<AGROW>
    iir(i).family = 'IIR';             %#ok<AGROW>
    iir(i).order  = 2*n;               %#ok<AGROW>
    iir(i).b      = [];                %#ok<AGROW>
    iir(i).sos    = sos;               %#ok<AGROW>
    iir(i).g      = g;                 %#ok<AGROW>
    fprintf('   %-16s max |pole| = %.5f  (stable: %d)\n', names{i}, ...
            max(abs(p{i})), all(abs(p{i}) < 1));
end

plot_filter_set(iir, fs, cfg, 'IIR band-pass filters (prototype order %d)', n, '06_iir_responses');

% pole-zero plots and group delay
fig = new_figure(cfg, 'IIR poles/zeros', 1200, 700);
for i = 1:3
    subplot(2, 3, i);
    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'k:'); hold on;
    plot(real(z{i}), imag(z{i}), 'bo', 'MarkerSize', 6);
    plot(real(p{i}), imag(p{i}), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold off; axis equal; grid on; xlim([-1.3 1.3]); ylim([-1.3 1.3]);
    title([names{i} ' pole-zero']); xlabel('Re'); ylabel('Im');

    subplot(2, 3, i+3);
    [H, f] = filter_response(iir(i), fs);
    ph = unwrap(angle(H));
    gd = -diff(ph) ./ diff(2*pi*f/fs);           % group delay in samples
    plot(f(2:end), gd, 'LineWidth', 1.2); grid on;
    xlim([0 fs/2]); ylim([0 max(5, min(200, 1.2*max(gd)))]);
    xlabel('Frequency (Hz)'); ylabel('Samples');
    title('Group delay (non-linear phase)');
end
save_figure(fig, cfg, '07_iir_polezero_groupdelay');

report_specs(iir, fs, cfg, 'task3_iir_specs.txt');

%% Analog prototype -> digital: bilinear transform vs impulse invariance
% (lab experiment: buttap / butter(...,'s'), bilinear, impinvar)
nd  = 4;
Om  = 2*fs*tan(pi*cfg.band/fs);               % pre-warped analog band edges (rad/s)
[bs, as]   = butter(nd, Om, 's');              % analog Butterworth band-pass
[bb, ab]   = bilinear_compat(bs, as, fs);      % bilinear, with pre-warping
[bs0, as0] = butter(nd, 2*pi*cfg.band, 's');   % analog design WITHOUT pre-warping
[bb0, ab0] = bilinear_compat(bs0, as0, fs);
[bi, ai]   = impinvar(bs0, as0, fs);           % impulse invariance
[bd, ad]   = butter(nd, cfg.band/(fs/2));      % built-in digital design

f  = linspace(1, fs/2, 4000)';
Ha = freqs(bs0, as0, 2*pi*f);                  % analog reference (true Hz axis)
Hb = freqz(bb, ab, f, fs);  Hb0 = freqz(bb0, ab0, f, fs);
Hi = freqz(bi, ai, f, fs);  Hd  = freqz(bd, ad, f, fs);
dB = @(H) 20*log10(abs(H(:)) + 1e-12);

fig = new_figure(cfg, 'Bilinear vs impulse invariance', 1150, 700);
subplot(2, 1, 1);
plot(f, dB(Ha), 'k', f, dB(Hb), 'b', f, dB(Hb0), 'm--', f, dB(Hi), 'r', 'LineWidth', 1.2);
grid on; xlim([0 fs/2]); ylim([-80 5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Analog prototype', 'Bilinear + pre-warping', 'Bilinear, no pre-warping', ...
       'Impulse invariance', 'Location', 'southwest');
title(sprintf('Butterworth band-pass (prototype order %d): analog-to-digital mapping', nd));
subplot(2, 1, 2);
plot(f, dB(Ha), 'k', f, dB(Hb), 'b', f, dB(Hb0), 'm--', f, dB(Hi), 'r', 'LineWidth', 1.2);
hold on; plot(cfg.band(2)*[1 1], [-20 2], 'k:'); hold off;
grid on; xlim([cfg.band(2)-1500 cfg.band(2)+2500]); ylim([-20 2]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Upper band edge (dotted line = specified 3800 Hz)');
save_figure(fig, cfg, '07b_iir_bilinear_vs_impinvar');

edge = @(H) f(find(f > mean(cfg.band) & dB(H) < -3, 1));
fid = fopen(fullfile(cfg.out_dir, 'task3_iir_specs.txt'), 'a');
l = sprintf(['\nAnalog-to-digital mapping (Butterworth, prototype order %d)\n' ...
             '  upper -3 dB edge, analog prototype      : %6.0f Hz\n' ...
             '  bilinear with pre-warping              : %6.0f Hz\n' ...
             '  bilinear without pre-warping           : %6.0f Hz (frequency warping)\n' ...
             '  impulse invariance                     : %6.0f Hz, attenuation at fs/2 = %.1f dB (aliasing)\n' ...
             '  max |coef diff| bilinear vs butter()   : %.2e (MATLAB designs digital IIR this way)\n'], ...
            nd, edge(Ha), edge(Hb), edge(Hb0), edge(Hi), -dB(Hi(end)), ...
            max([abs(bb(:) - bd(:)); abs(ab(:) - ad(:))]));
fprintf('%s', l); fprintf(fid, '%s', l);
fclose(fid);
end
