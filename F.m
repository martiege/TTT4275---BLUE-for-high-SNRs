function res = F(N, w_0_est, T)
res = 0;

for n = 1:N
    res = res + exp(1j * w_0_est * n * T);
end

res = res / N;
end

