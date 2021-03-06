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
var = (A^2 / 2) ./ db2mag(SNR);
est = zeros(length(SNR), 2);
H_b = [T * n', ones(length(n), 1)];
C_base = eye(N);


for i = 1:length(SNR)
    % variance = A^2 / (2 * db2mag(SNR(i)));
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(var(i)));
    est(i, :) = BLUE(x(i, :), H_b, var(i) * C_base);
    % BLUE(x(i, :), H, var(i) * C_base);
end

CRLB_omega = (12 / (A^2 * T^2 * N * (N^2 - 1))) .* var;
CRLB_phi = (12 * (n_0^2 * N + 2 * n_0 * P + Q) / (A^2 * N^2 * (N^2 - 1))) .* var;
BLUE_omega = zeros(1, length(SNR));
BLUE_phi = zeros(1, length(SNR));
for i = 1:length(SNR)
    C = inv(H_b' * inv(var(i) * C_base) * H_b);
    BLUE_omega(1, i) = C(1, 1);
    BLUE_phi(1, i) = C(2, 2);
end

differ = zeros(length(SNR), N);
for snr = 1:length(SNR)
    for i = 1:N-1
        differ(snr, i) = angle(x(snr, i + 1)) - angle(x(snr, i));
    end
end

SNR = -10:10:40;
differ = zeros(length(SNR), N - 1);
var = (A^2 / 2) ./ db2mag(SNR);
est2 = zeros(length(SNR), 2);
H_c_base = T * ones(N-1, 1);
D_base = diag([ones(1, N-1), 0]);
D_base = D_base(1:N-1, 1:N);
D = circshift(D_base, 1, 2) - D_base;
H_c = H_c_base;
C_c_base = D * C_base * D';

for i = 1:length(SNR)
    y = angle(x(i, :));
    for j = 1:N-1
        differ(i, j) = mod(y(j + 1) - y(j), pi);
    end
    
    C_c = var(i) * C_c_base;
    est2(i, 1) = inv(H_c' * inv(C_c) * H_c) * H_c' * inv(C_c) * differ(i, :)';
    F_sum = 0;
    for j = 1:N
        F_sum = F_sum + x(i, j) * exp(-1j * est2(i, 1) * n(i) * T);
    end
    
    est2(i, 2) = mod(angle(exp(-1j * est2(i, 1) * n_0 * T) * (F_sum / N)), pi);
end




