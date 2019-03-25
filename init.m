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
% x = A * exp(j * (w_0 * n * T + phi));
% v =  -1 + 2 .* rand(N,s1) + j(- 1 + 2 .* rand(N,1));
%SNR = A^2 / (2 * sigma^2);

SNR = -10:10:40;
x = zeros(length(SNR), N);
for i = 1:length(SNR)
    variance = A^2 / (2 * db2mag(SNR(i)));
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(variance));
end


H = ones(N, 2);
for i = 1:N 
    H(i, 1) = i * T;
end

C = A^2 / (2 * db2mag(SNR(length(SNR)))) * eye(N);
est = (inv(H' * inv(C) * H) * H' * inv(C)) * unwrap(angle(x(length(SNR), :)))';

w_0_est = est(1);
phi_est = est(2);

x_gen = unwrap(angle(x(length(SNR), :)));
x_est = unwrap(angle(A * exp(1j * (w_0_est * n' * T + phi_est))));

plot(n, x_gen, n, x_est);
legend("x_gen", "x_est");

% plot(n, unwrap(angle(x(1, :))), n, unwrap(angle(x(2, :))), n, unwrap(angle(x(3, :))), n, unwrap(angle(x(4, :))), n, unwrap(angle(x(4, :))), n, unwrap(angle(x(6, :))));
% legend("SNR -10", "SNR 0" , "SNR 10" , "SNR 20" , "SNR 30", "SNR 40");