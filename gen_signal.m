function x = gen_signal(omega_0, n, A, T, phi, w_e, w_sigma)
%GEN_ Summary of this function goes here
%   Detailed explanation goes here

% w = (rand(N, 1) + 1j * rand(N, 1)) .* (2 * w_sigma) - (w_sigma - w_e);
w = normrnd(w_e, w_sigma) + 1j * normrnd(w_e, w_sigma);

x = A * exp(1j * (omega_0 * n' * T + phi)) + w;

end

