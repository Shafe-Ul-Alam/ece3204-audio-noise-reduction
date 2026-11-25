function [bz, az] = bilinear_compat(bs, as, fs)
% BILINEAR_COMPAT  Bilinear transform s = 2 fs (z-1)/(z+1) that runs in MATLAB and Octave.
%   MATLAB : bilinear(b, a, fs)   (third argument = sampling frequency)
%   Octave : bilinear(b, a, T)    (third argument = sampling period)
if exist('OCTAVE_VERSION', 'builtin')
    [bz, az] = bilinear(bs, as, 1/fs);
else
    [bz, az] = bilinear(bs, as, fs);
end
bz = bz / az(1); az = az / az(1);
end
