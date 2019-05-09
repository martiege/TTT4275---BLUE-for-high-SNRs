clc; clear all; close all;

%% Declaring variables
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
var = (A^2 / 2) ./ db2mag(SNR);

CRLB_omega = (12 / (A^2 * T^2 * N * (N^2 - 1))) .* var;
CRLB_phi = (12 * (n_0^2 * N + 2 * n_0 * P + Q) / ... 
            (A^2 * N^2 * (N^2 - 1))) .* var;


x = zeros(length(SNR), N);
v = zeros(length(SNR), N);
x_diff = zeros(length(SNR), N - 1);


% for b)
est_b = zeros(length(SNR), 2);
H_b = [T * n', ones(length(n), 1)];
C_base = eye(N);

%for c)
est_c = zeros(length(SNR), 2);
H_c = T * ones(N-1, 1);
D_base = diag([ones(1, N-1), 0]);
D_base = D_base(1:N-1, 1:N);
D = circshift(D_base, 1, 2) - D_base;

%% Estimating per SNR
for i = 1:length(SNR)
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, 0);
    v(i, :) = gen_noise(N, 0, sqrt(var(i)));
    
    y = angle(x(i, :))';
    for j = 1:N-1
        x_diff(i, j) = w_0*T + v(j+1) - v(j);
    end
    
    C = D * (var(i)) * D';
    
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(var(i)));
    
    % Estimates in b) (phi is unwrapped)
    est_b(i, :) = BLUE(x(i, :), H_b, var(i) * C_base);
    
    % Omega estimate in c)
    est_c(i, 1) = abs(BLUE_c(x_diff(i, :)', H_c, C));
    
    
    % MLE for phase estimate in c)
    res = zeros(N, 1);
    for n_ = 1:N
        res(n_) = x(i, n_) * exp(-1j * est_c(i, 1) * (n_ - 1)*T);
    end
    
    fourier = mean(res);
    est_c(i, 2) = angle(exp(-1j*est_c(i,1)*n_0*T)*fourier);
    
end