function [x, fs] = load_audio(cfg)
% LOAD_AUDIO  Read the reference WAV, convert to mono, trim, resample, normalise.
%
%   [x, fs] = load_audio(cfg) returns a column vector x at cfg.fs Hz.

file = cfg.input_file;
if ~exist(file, 'file')
    % fall back to the first audio file found in the audio/ folder
    [folder, ~, ~] = fileparts(file);
    cand = [dir(fullfile(folder, '*.wav')); dir(fullfile(folder, '*.flac')); ...
            dir(fullfile(folder, '*.mp3')); dir(fullfile(folder, '*.ogg'))];
    if isempty(cand)
        error(['No input audio found. Put a clean speech recording in the audio/ ' ...
               'folder, e.g. audio/clean_speech.wav (see README.md).']);
    end
    file = fullfile(folder, cand(1).name);
end

[x, fs_in] = audioread(file);

% 1) stereo -> mono (average channels)
if size(x, 2) > 1
    x = mean(x, 2);
end

% 2) keep at most cfg.max_duration seconds
n_max = round(cfg.max_duration * fs_in);
if length(x) > n_max
    x = x(1:n_max);
end

% 3) resample to the working rate (rational factor p/q)
if fs_in ~= cfg.fs
    [p, q] = rat(cfg.fs / fs_in);
    x = resample(x, p, q);
end
fs = cfg.fs;

% 4) remove DC offset and normalise peak
x = x - mean(x);
x = cfg.peak_level * x / max(abs(x));

fprintf('Loaded %s : %.2f s, original fs = %d Hz -> working fs = %d Hz\n', ...
        file, length(x)/fs, fs_in, fs);
end
