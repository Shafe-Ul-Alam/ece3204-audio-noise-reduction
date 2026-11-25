function h = fir_window_design(N, fc, fs, ftype, win)
% FIR_WINDOW_DESIGN  User-defined FIR design by the window method (Li Tan, Ch. 7).
%
%   h = fir_window_design(N, fc, fs, ftype, win)
%     N     : number of taps, odd (N = 2M + 1)
%     fc    : cut-off frequency in Hz (scalar) or [fl fh] for band filters
%     fs    : sampling frequency (Hz)
%     ftype : 'lowpass' | 'highpass' | 'bandpass' | 'bandstop'
%     win   : 'rectangular' | 'hanning' | 'hamming' | 'blackman'
%
%   Step 1 - ideal (truncated) impulse response, k = -M..M, h(-k) = h(k):
%     lowpass  : h(0) = wc/pi,           h(k) = sin(k wc)/(k pi)
%     highpass : h(0) = (pi - wc)/pi,    h(k) = -sin(k wc)/(k pi)
%     bandpass : h(0) = (wh - wl)/pi,    h(k) = sin(k wh)/(k pi) - sin(k wl)/(k pi)
%     bandstop : h(0) = (pi - wh + wl)/pi, h(k) = -sin(k wh)/(k pi) + sin(k wl)/(k pi)
%   Step 2 - multiply by the window w(k), k = -M..M:
%     rectangular : 1
%     hanning     : 0.5  + 0.5 cos(k pi / M)
%     hamming     : 0.54 + 0.46 cos(k pi / M)
%     blackman    : 0.42 + 0.5 cos(k pi / M) + 0.08 cos(2 k pi / M)
%   The result is causal after a delay of M samples and has linear phase.

if mod(N, 2) == 0
    N = N + 1;                       % force odd length (type-I linear phase)
end
M  = (N - 1) / 2;
k  = (-M:M)';
w  = 2*pi*fc(:)'/fs;                 % digital cut-off(s) in rad/sample

% ---- step 1: ideal impulse response ----
hd = zeros(N, 1);
nz = k ~= 0;
switch lower(ftype)
    case 'lowpass'
        hd(k == 0) = w(1)/pi;
        hd(nz) = sin(k(nz)*w(1)) ./ (k(nz)*pi);
    case 'highpass'
        hd(k == 0) = (pi - w(1))/pi;
        hd(nz) = -sin(k(nz)*w(1)) ./ (k(nz)*pi);
    case 'bandpass'
        hd(k == 0) = (w(2) - w(1))/pi;
        hd(nz) = sin(k(nz)*w(2)) ./ (k(nz)*pi) - sin(k(nz)*w(1)) ./ (k(nz)*pi);
    case 'bandstop'
        hd(k == 0) = (pi - w(2) + w(1))/pi;
        hd(nz) = -sin(k(nz)*w(2)) ./ (k(nz)*pi) + sin(k(nz)*w(1)) ./ (k(nz)*pi);
    otherwise
        error('Unknown filter type: %s', ftype);
end

% ---- step 2: window ----
switch lower(win)
    case {'rectangular', 'rect'}
        wn = ones(N, 1);
    case {'hanning', 'hann'}
        wn = 0.5 + 0.5*cos(k*pi/M);
    case 'hamming'
        wn = 0.54 + 0.46*cos(k*pi/M);
    case 'blackman'
        wn = 0.42 + 0.5*cos(k*pi/M) + 0.08*cos(2*k*pi/M);
    otherwise
        error('Unknown window: %s', win);
end

h = (hd .* wn)';                     % row vector, like fir1
end
