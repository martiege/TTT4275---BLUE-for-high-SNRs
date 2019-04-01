clc; clear all; close all;

F_s = 10^6;
T = 1 / F_s;
f_0 = 10^5;
w_0 = 2 * pi * f_0;
phi = pi / 8;
A = 1;
N = 513;
P = N * (N - 1) / 2;
Q = N * (N - 1) * (2 * N - 1) / 6;
n_0 = -P / N;
n_N = n_0 + N - 1;
n = n_0:n_N;

SNR = -10:10:40;
x = zeros(length(SNR), N);
diff = zeros(length(SNR), N - 1);
var = (A^2 / 2) ./ db2mag(SNR);
est = zeros(length(SNR), 2);
H = T * ones(N-1, 1);
D_base = diag([ones(1, N-1), 0]);
D_base = D_base(1:N-1, 1:N);
D = circshift(D_base, 1, 2) - D_base;

for i = 1:length(SNR)
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(var(i)));
    for j = 1:N-1
        diff(i, j) = angle(x(i, j + 1)) - angle(x(i, j));
    end
    
    C = D * (var(i)) * D';
    est(i, 1) = BLUE_c(diff(i, :)', H, C);
    est(i, 1) = BLUE_c(D * angle(x(i, :)'), H, C);
    % C_inv = inv(var(i) * (D * D'));
    % est(i, 1) = inv(H' * C_inv * H) * (H' * C_inv * (D * angle(x(i, :)')));
    %inv(H' * inv(var(i) * D * D') * H) * H' * inv(var(i) * D * D') * D * angle(x(i, :)');
    %BLUE_c(D * angle(x(i, :)'), H, D * var(i) * eye(N) * D');
    %est(i, 2) = (1 / N) * sum(angle(x(i, :) - est(i, 1) * n));
end



