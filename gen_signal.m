function x = gen_signal(omega_0, N, A, T, phi, w_e, w_sigma)
%GEN_ Summary of this function goes here
%   Detailed explanation goes here

w = (rand(N, 1) + 1j * rand(N, 1)) .* (2 * w_sigma) - (1 - w_e);

N = 513;
P = N * (N - 1) / 2;
Q = N * (N - 1) * (2 * N - 1) / 6;
n_0 = -P / N;
n_N = n_0 + N - 1;
n = n_0:n_N;

x = A * exp(1j * (omega_0 * n' * T + phi)) + w;

end

