function s = filter_specs(flt, fs, cfg)
% FILTER_SPECS  Measure the achieved specifications of a designed filter
% (single pass; forward-backward filtering doubles every dB value).
%   s.ripple_db : peak-to-peak pass-band deviation (dB) in the inner pass-band 300 Hz - (f2 - 300 Hz)
%   s.atten_db  : minimum attenuation (dB) in the upper stop-band f >= f2 + 700 Hz
%   s.hum_db    : attenuation at the 50 Hz mains frequency (dB)
%   s.trans_hz  : upper transition width, from the -1 dB to the -40 dB point (Hz)
%   s.n_mult    : multiplications per output sample (computational cost)
[H, f] = filter_response(flt, fs);
M = 20*log10(abs(H) + 1e-12);
pb = f >= 300 & f <= cfg.band(2) - 300;
sb = f >= cfg.band(2) + 700;
s.ripple_db = max(M(pb)) - min(M(pb));
s.atten_db  = -max(M(sb));
[~, i50] = min(abs(f - cfg.hum_f0));
s.hum_db    = -M(i50);
up  = f > cfg.band(2) - 500;
f1  = f(find(up & M < -1, 1));
f40 = f(find(up & M < -40, 1));
if isempty(f40) || isempty(f1), s.trans_hz = NaN; else, s.trans_hz = f40 - f1; end
switch flt.type
    case 'fir', s.n_mult = length(flt.b);
    case 'iir', s.n_mult = 5 * size(flt.sos, 1);
end
end
