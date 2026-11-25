function noisy = make_noisy_set(clean, fs, cfg)
% MAKE_NOISY_SET  Create several noisy versions of the clean signal.
%
%   noisy(k).name, noisy(k).x, noisy(k).noise, noisy(k).snr_target
%
%   Noise types (each scaled to the target input SNR in cfg.noise_snr):
%     white : wide-band Gaussian noise (e.g. microphone / thermal noise)
%     hum   : 50 Hz mains hum + harmonics (electrical interference)
%     hiss  : high-frequency noise above cfg.hiss_cutoff + a 6 kHz tone
%     mixed : hum + hiss + white together (realistic worst case)

rng(cfg.seed);                       % reproducible noise
N = length(clean);
t = (0:N-1)' / fs;

% ---- building blocks (un-scaled) ----
white = randn(N, 1);

hum = zeros(N, 1);
for k = 1:cfg.hum_harm                       % harmonics decay as 1/k
    hum = hum + (1/k) * sin(2*pi*k*cfg.hum_f0*t + 2*pi*rand);
end

b_hp = fir_window_design(129, cfg.hiss_cutoff, fs, 'highpass', 'hamming'); % shape white noise into hiss
hiss = filter(b_hp, 1, randn(N, 1));
hiss = hiss / rms(hiss) + 0.7 * sqrt(2) * sin(2*pi*cfg.tone_freq*t);

blocks.white = white;
blocks.hum   = hum;
blocks.hiss  = hiss;
blocks.mixed = hum/rms(hum) + hiss/rms(hiss) + 0.5*white;

% ---- scale each noise to the required SNR ----
P_s = mean(clean.^2);
for k = 1:numel(cfg.noise_types)
    name = cfg.noise_types{k};
    n    = blocks.(name);
    snr_db = cfg.noise_snr(k);
    n = n * sqrt(P_s / (mean(n.^2) * 10^(snr_db/10)));

    noisy(k).name       = name;            %#ok<AGROW>
    noisy(k).noise      = n;               %#ok<AGROW>
    noisy(k).x          = clean + n;       %#ok<AGROW>
    noisy(k).snr_target = snr_db;          %#ok<AGROW>
end
end

function r = rms(v)
r = sqrt(mean(v.^2));
end
