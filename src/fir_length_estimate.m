function N = fir_length_estimate(win, df, fs)
% FIR_LENGTH_ESTIMATE  Required number of taps for a window design (Li Tan, Table 7.7).
%
%   N = factor * fs / df, rounded up to the next odd integer
%     window       factor   pass-band ripple   stop-band attenuation
%     rectangular   0.9        0.7416 dB             21 dB
%     hanning       3.1        0.0546 dB             44 dB
%     hamming       3.3        0.0194 dB             53 dB
%     blackman      5.5        0.0017 dB             74 dB
%   df = transition width in Hz (|stop-band edge - pass-band edge|).

switch lower(win)
    case {'rectangular', 'rect'}, c = 0.9;
    case {'hanning', 'hann'},     c = 3.1;
    case 'hamming',               c = 3.3;
    case 'blackman',              c = 5.5;
    otherwise, error('Unknown window: %s', win);
end
N = ceil(c * fs / df);
if mod(N, 2) == 0
    N = N + 1;
end
end
