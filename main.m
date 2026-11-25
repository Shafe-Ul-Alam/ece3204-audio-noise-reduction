%% =======================================================================
%  ECE 3204 - Open-Ended DSP Laboratory Project
%  Design and Performance Evaluation of an Adaptive Audio Noise Reduction
%  System Using Digital Signal Processing Techniques
%
%  Run this file from the project root folder. It executes all 7 tasks:
%    1. WAV loading, waveform and FFT analysis         (DFT / FFT)
%    2. FIR filters: Rectangular / Hamming / Blackman  (FIR design)
%    3. IIR filters: Butterworth / Chebyshev / Elliptic(IIR design)
%    4. Down/up-sampling, aliasing and imaging         (Multirate DSP)
%    5. Periodogram vs Welch PSD                        (PSD)
%    6. DWT wavelet denoising                          (Wavelet)
%    7. Evaluation, comparison and parameter sweeps    (Integrated)
%
%  All parameters live in src/config.m. Results go to results/.
%  =======================================================================
clear; close all; clc;
addpath('src');

% Octave compatibility (MATLAB ignores this block)
if exist('OCTAVE_VERSION', 'builtin')
    pkg load signal;
end

cfg = config();
for d = {cfg.out_dir, cfg.fig_dir, cfg.audio_dir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

t_start = tic;

%% Input and noisy versions
[clean, fs] = load_audio(cfg);
noisy = make_noisy_set(clean, fs, cfg);

%% Tasks
task1_analysis(clean, noisy, fs, cfg);
fir = task2_fir_design(fs, cfg);
iir = task3_iir_design(fs, cfg);
task4_multirate(clean, fs, cfg);
task5_psd(clean, noisy, fs, cfg);
task6_wavelet(clean, noisy, fs, cfg);
results = task7_evaluate(clean, noisy, fs, [fir, iir], cfg);
task7_sweeps(clean, noisy, fs, cfg);

fprintf('\nAll tasks finished in %.1f s. See the results/ folder.\n', toc(t_start));
