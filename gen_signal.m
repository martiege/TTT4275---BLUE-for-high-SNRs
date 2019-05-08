function x = gen_signal(omega_0, n, A, T, phi, w_e, w_sigma)

w = gen_noise(length(n), w_e, w_sigma);
x = A * exp(1j * (omega_0 * n * T + phi)) + w;

end

