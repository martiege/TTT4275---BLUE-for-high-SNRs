function v = gen_noise(N, w_e, w_sigma)
v = normrnd(w_e, w_sigma, 1, N) + 1j * normrnd(w_e, w_sigma, 1, N);

end