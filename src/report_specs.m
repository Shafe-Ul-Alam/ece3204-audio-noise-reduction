function report_specs(flts, fs, cfg, fname)
% REPORT_SPECS  Print and save the measured specifications of a filter set.
fid = fopen(fullfile(cfg.out_dir, fname), 'w');
hdr = sprintf('%-18s %6s %11s %12s %11s %15s %10s\n', 'Filter', 'Order', ...
              'PassDev(dB)', 'StopAtt(dB)', '@50Hz(dB)', 'Transition(Hz)', 'Mult/samp');
fprintf('%s', hdr); fprintf(fid, '%s', hdr);
for k = 1:numel(flts)
    s = filter_specs(flts(k), fs, cfg);
    line = sprintf('%-18s %6d %11.3f %12.1f %11.1f %15.0f %10d\n', flts(k).name, ...
                   flts(k).order, s.ripple_db, s.atten_db, s.hum_db, s.trans_hz, s.n_mult);
    fprintf('%s', line); fprintf(fid, '%s', line);
end
fprintf(fid, ['\nValues are for a single pass. The system filters forward and backward\n' ...
              '(zero phase), so effective attenuation in dB is doubled.\n']);
fclose(fid);
end
