function [H, f] = filter_response(flt, fs, nfreq)
% FILTER_RESPONSE  Complex frequency response of a designed FIR/IIR filter.
%   For IIR filters the response is the product of all second-order sections.
if nargin < 3, nfreq = 8192; end
switch flt.type
    case 'fir'
        [H, f] = freqz(flt.b, 1, nfreq, fs);
    case 'iir'
        H = flt.g * ones(nfreq, 1);
        for s = 1:size(flt.sos, 1)
            [Hs, f] = freqz(flt.sos(s, 1:3), flt.sos(s, 4:6), nfreq, fs);
            H = H .* Hs;
        end
end
H = H(:); f = f(:);
end
