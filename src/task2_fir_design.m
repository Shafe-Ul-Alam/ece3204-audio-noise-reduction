function fir = task2_fir_design(fs, cfg)
% TASK 2 - FIR band-pass filters by the window method (user-defined design).
%
% The filters are designed with fir_window_design.m, which implements the
% Li Tan window method used in the lab (ideal impulse response x window),
% for the Rectangular, Hanning, Hamming and Blackman windows. The window
% trades main-lobe width (transition band) against side-lobe level
% (stop-band attenuation):
%   Rectangular : narrowest transition, ~21 dB stop-band, Gibbs ripple
%   Hanning     : ~44 dB stop-band
%   Hamming     : ~53 dB stop-band, moderate transition width
%   Blackman    : ~74 dB stop-band, widest transition
% All designs have exactly linear phase (symmetric coefficients).

fprintf('\n=== Task 2: FIR filter design (user-defined window method) ===\n');
Ntaps = cfg.fir_order + 1;
wins  = cfg.fir_windows;
fid = fopen(fullfile(cfg.out_dir, 'task2_fir_design.txt'), 'w');

% ---- required length from the lab formula N = c*fs/df ----
l = sprintf(['Required taps from N = c*fs/df (Li Tan) for transition width df = %d Hz:\n'], cfg.fir_df);
fprintf('%s', l); fprintf(fid, '%s', l);
for k = 1:numel(wins)
    l = sprintf('   %-12s N = %d\n', wins{k}, fir_length_estimate(wins{k}, cfg.fir_df, fs));
    fprintf('%s', l); fprintf(fid, '%s', l);
end
l = sprintf(['Used: N = %d taps for every window (equal cost, fair comparison).\n' ...
             'The lower edge (%d Hz) is close to 0 Hz, so the lower transition band\n' ...
             'is the hardest constraint for the wide-main-lobe windows.\n\n'], Ntaps, cfg.band(1));
fprintf('%s', l); fprintf(fid, '%s', l);

for k = 1:numel(wins)
    fir(k).name   = ['FIR-' upper(wins{k}(1)) wins{k}(2:end)];          %#ok<AGROW>
    fir(k).type   = 'fir';                                               %#ok<AGROW>
    fir(k).family = 'FIR';                                               %#ok<AGROW>
    fir(k).order  = Ntaps - 1;                                           %#ok<AGROW>
    fir(k).b      = fir_window_design(Ntaps, cfg.band, fs, 'bandpass', wins{k}); %#ok<AGROW>
    fir(k).sos    = [];                                                  %#ok<AGROW>
    fir(k).g      = 1;                                                   %#ok<AGROW>
end

% ---- cross-check with MATLAB's built-in fir1 (Hamming) ----
try
    b_ref = fir1(Ntaps - 1, cfg.band/(fs/2), hamming(Ntaps));
    iH = find(strcmp(wins, 'hamming'));
    [H1, ff] = freqz(fir(iH).b, 1, 4096, fs); H2 = freqz(b_ref, 1, 4096, fs);
    pb = ff >= 300 & ff <= cfg.band(2) - 300;
    l = sprintf(['Check vs built-in fir1 (Hamming): max pass-band difference = %.4f dB,\n' ...
                 'max coefficient difference = %.2e\n' ...
                 '(fir1 additionally scales the coefficients for exactly 0 dB at the pass-band centre)\n\n'], ...
                max(abs(20*log10(abs(H1(pb))) - 20*log10(abs(H2(pb))))), max(abs(fir(iH).b - b_ref)));
    fprintf('%s', l); fprintf(fid, '%s', l);
catch
    fprintf('   (fir1 cross-check skipped)\n');
end
fclose(fid);

plot_filter_set(fir, fs, cfg, 'FIR band-pass filters (window method, %d taps)', Ntaps, '04_fir_responses');

% impulse responses
fig = new_figure(cfg, 'FIR impulse responses', 1200, 420);
for k = 1:numel(fir)
    subplot(1, numel(fir), k);
    stem(0:Ntaps-1, fir(k).b, 'Marker', 'none'); grid on;
    xlim([0 Ntaps-1]); xlabel('n'); ylabel('h[n]');
    title(fir(k).name);
end
save_figure(fig, cfg, '05_fir_impulse_responses');

report_specs(fir, fs, cfg, 'task2_fir_specs.txt');
end
