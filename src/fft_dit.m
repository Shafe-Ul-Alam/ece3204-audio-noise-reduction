function X = fft_dit(x)
% FFT_DIT  Radix-2 decimation-in-time FFT (recursive), N must be a power of 2.
%   Split x into even and odd samples (each an N/2-point DFT) and combine
%   with the butterfly
%       X(k)       = E(k) + W_N^k O(k)
%       X(k + N/2) = E(k) - W_N^k O(k),   k = 0 .. N/2-1
%   Cost: (N/2) log2(N) complex multiplications, O(N log N).
x = x(:);
N = length(x);
if N == 1
    X = x;
    return;
end
if mod(N, 2) ~= 0
    error('fft_dit: length must be a power of 2');
end
E = fft_dit(x(1:2:end));                 % even-indexed samples
O = fft_dit(x(2:2:end));                 % odd-indexed samples
W = exp(-1j*2*pi*(0:N/2-1)'/N);          % twiddle factors
X = [E + W.*O; E - W.*O];
end
