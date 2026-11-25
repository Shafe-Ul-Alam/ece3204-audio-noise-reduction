function write_wav(y, fs, cfg, fname)
% WRITE_WAV  Save audio to results/audio, scaling down only if it would clip.
%   audiowrite needs an integer sample rate. After decimation by 3 the rate
%   is 16000/3 = 5333.33 Hz, so it is rounded (5333 Hz, 0.006 % pitch error,
%   inaudible). Processing itself always uses the exact rate.
y = y(:);
pk = max(abs(y));
if pk > 0.99
    y = 0.99 * y / pk;
end
audiowrite(fullfile(cfg.audio_dir, [fname '.wav']), y, round(fs));
end
