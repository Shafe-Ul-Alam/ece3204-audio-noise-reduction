function X = dft_direct(x)
% DFT_DIRECT  DFT computed straight from its definition (no fft()).
%   X(k) = sum_{n=0}^{N-1} x(n) W_N^{kn},  W_N = exp(-j 2 pi / N)
%   Written as a matrix-vector product: N^2 complex multiplications, O(N^2).
x = x(:);
N = length(x);
n = 0:N-1;
W = exp(-1j*2*pi*(n')*n/N);      % N x N twiddle-factor matrix W_N^{kn}
X = W * x;
end
