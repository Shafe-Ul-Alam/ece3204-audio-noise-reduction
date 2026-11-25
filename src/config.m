function cfg = config()
% CONFIG  All project parameters in one place (reproducible experiments).
%
%   cfg = config() returns a struct used by every task script.
%   Change values here, not inside the task files.

%% ---------------- Input audio ----------------
cfg.input_file   = fullfile('audio', 'clean_speech.wav'); % clean reference clip
cfg.fs           = 16000;   % working sample rate (Hz) - standard for speech
cfg.max_duration = 10;      % seconds of audio used (keeps runtime short)
cfg.peak_level   = 0.5;     % clean signal is normalised to this peak

%% ---------------- Noise generation ----------------
cfg.seed        = 42;                               % reproducible noise
cfg.noise_types = {'white', 'hum', 'hiss', 'mixed'};
cfg.noise_snr   = [5, 0, 0, 0];                     % input SNR (dB) per type
cfg.hum_f0      = 50;                               % mains frequency in Bangladesh (Hz)
cfg.hum_harm    = 5;                                % number of hum harmonics
cfg.hiss_cutoff = 4500;                             % hiss = noise above this (Hz)
cfg.tone_freq   = 6000;                             % high-frequency interference tone (Hz)

%% ---------------- Filter specifications ----------------
% Speech energy that matters for intelligibility lies roughly in 150-3800 Hz.
% A band-pass filter therefore removes hum (<150 Hz) and hiss/tone (>3.8 kHz).
cfg.band       = [150 3800];   % pass-band edges (Hz)
cfg.fir_order  = 256;          % FIR order (taps = order + 1)
cfg.fir_windows = {'rectangular', 'hanning', 'hamming', 'blackman'};
cfg.fir_df     = 150;          % transition width used for the length estimate (Hz)
cfg.iir_order  = 6;            % IIR prototype order (band-pass => 2*order)
cfg.cheby_Rp   = 1;            % Chebyshev-I / Elliptic pass-band ripple (dB)
cfg.ellip_Rs   = 60;           % Elliptic stop-band attenuation (dB)

%% ---------------- Parameter sweeps (design challenge 2) ----------------
cfg.sweep_fir_orders = [32 64 128 256 512];
cfg.sweep_iir_orders = [2 4 6 8 10];
cfg.sweep_upper_fc   = [2500 3000 3400 3800 4500 5500];

%% ---------------- Multirate ----------------
cfg.mr_factors   = [2 3 4];   % down/up-sampling factors to investigate
cfg.mr_fir_order = 128;       % anti-aliasing / interpolation filter order

%% ---------------- PSD ----------------
cfg.nfft          = 1024;
cfg.welch_win     = 512;      % segment length (samples)
cfg.welch_overlap = 0.5;      % 50 % overlap
cfg.ar_orders     = [12 40];   % AR model orders (Yule-Walker / Burg)

%% ---------------- Wavelet ----------------
cfg.wavelet       = 'db4';      % haar | db2 | db4 | db8 | sym4 | coif1
cfg.wav_families  = {'haar', 'db2', 'db4', 'db8', 'sym4', 'coif1'};
cfg.wav_levels    = 5;          % DWT decomposition levels
cfg.wav_threshold = 'soft';     % 'soft' or 'hard'
cfg.wav_rule      = 'bayes';    % 'bayes'  = BayesShrink (adaptive per level)
                                % 'global' = VisuShrink universal threshold
                                % 'level'  = universal threshold, sigma per level

%% ---------------- Adaptive (NLMS) hum canceller ----------------
cfg.nlms_mu    = 0.01;   % normalised step size
cfg.nlms_harm  = 5;      % harmonics in the reference signal

%% ---------------- Evaluation / output ----------------
cfg.seg_len       = 256;     % frame length for segmental SNR / LSD
cfg.timing_runs   = 3;       % runtime averaged over this many runs
cfg.show_figures  = false;   % true = figures pop up; false = only saved
cfg.out_dir       = 'results';
cfg.fig_dir       = fullfile('results', 'figures');
cfg.audio_dir     = fullfile('results', 'audio');
end
