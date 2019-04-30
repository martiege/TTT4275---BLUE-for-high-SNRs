function res = F(x, N, w_0_est, T)
%F Summary of this function goes here
%   Detailed explanation goes here
res = 0;

for n = 1:N
    res = res + x(n) * exp(-1j * w_0_est * n * T);
end

res = res / N;
end

