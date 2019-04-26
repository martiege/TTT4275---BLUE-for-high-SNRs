function v = gen_noise(N, w_e, w_sigma)
%v = (rand(1, N) + 1j * rand(1, N)) .* (2 * w_sigma) - (w_sigma - w_e);
v = normrnd(w_e, w_sigma, 1, N) + 1j * normrnd(w_e, w_sigma, 1, N);

end