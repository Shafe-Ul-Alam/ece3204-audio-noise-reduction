function y = apply_filter(flt, x)
% APPLY_FILTER  Zero-phase (forward-backward) filtering of x with a designed filter.
%
%   flt.type = 'fir' -> uses flt.b
%   flt.type = 'iir' -> uses flt.sos and flt.g (second-order sections, numerically
%                       stable for high-order band-pass designs)
%
%   Forward-backward filtering removes phase distortion and group delay, so the
%   output is time-aligned with the clean reference and SNR can be compared fairly
%   between FIR (linear phase) and IIR (non-linear phase) designs.
%   Note: the effective magnitude response becomes |H|^2.

x = x(:);
L = min(3000, length(x) - 1);               % reflection padding reduces edge transients
xp = [2*x(1) - flipud(x(2:L+1)); x; 2*x(end) - flipud(x(end-L:end-1))];

switch flt.type
    case 'fir'
        y = filter(flt.b, 1, xp);
        y = flipud(filter(flt.b, 1, flipud(y)));
    case 'iir'
        y = flt.g * sosfilt(flt.sos, xp);
        y = flipud(flt.g * sosfilt(flt.sos, flipud(y)));
    otherwise
        error('Unknown filter type: %s', flt.type);
end

y = y(L+1 : L+length(x));
end
