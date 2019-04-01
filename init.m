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
H = [T * n', ones(length(n), 1)];
C_base = eye(N);


for i = 1:length(SNR)
    % variance = A^2 / (2 * db2mag(SNR(i)));
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(var(i)));
    est(i, :) = BLUE(x(i, :), H, var(i) * C_base);
    % BLUE(x(i, :), H, var(i) * C_base);
end

CRLB_omega = (12 / (A^2 * T^2 * N * (N^2 - 1))) .* var;
CRLB_phi = (12 * (n_0^2 * N + 2 * n_0 * P + Q) / (A^2 * N^2 * (N^2 - 1))) .* var;
BLUE_omega = zeros(1, length(SNR));
BLUE_phi = zeros(1, length(SNR));
for i = 1:length(SNR)
    C = inv(H' * inv(var(i) * C_base) * H);
    BLUE_omega(1, i) = C(1, 1);
    BLUE_phi(1, i) = C(2, 2);
end

f1 = figure(1);
grid on;
title('Variance of estimated omega for CRLB and BLUE');
xlabel('SNR [dB]');
ylabel('Variance of the estimate of omega from CRLB and BLUE');
plot(SNR, CRLB_omega, SNR, BLUE_omega);
legend('CRLB omega', 'BLUE omega');

saveas(f1, 'figures/var_est_omega', 'epsc');

f2 = figure(2);
grid on;
title('Difference between the variance of the estimated omega');
xlabel('SNR [dB]');
ylabel('Difference of the variances of the estimate of omega from CRLB and BLUE');
plot(SNR, BLUE_omega - CRLB_omega);
legend('Omega: Difference between CRLB and BLUE');

saveas(f2, 'figures/diff_var_est_omega', 'epsc');

f3 = figure(3);
grid on;
title('Plot of the variance of the estimated phi given by CRLB and BLUE');
xlabel('SNR [dB]');
ylabel('Variance of the estimate of phi from CRLB and BLUE');
plot(SNR, CRLB_phi, SNR, BLUE_phi);
legend('CRLB phi', 'BLUE phi');

saveas(f3, 'figures/var_est_phi', 'epsc');

f4 = figure(4);
grid on;
title('Difference between the variance of the estimated phi');
xlabel('SNR [dB]');
ylabel('Difference of the variances of the estimate of phi from CRLB and BLUE');
plot(SNR, BLUE_phi - CRLB_phi);
legend('Phi: Difference between CRLB and BLUE');

saveas(f4, 'figures/diff_var_est_phi', 'epsc');

figure(5);
grid on;
for i = 1:length(SNR)
   sig = A * exp(1j * (est(i, 1) * n' * T + mod(est(i, 2), pi)));
   y = unwrap(angle(sig));
   plot(n, y)
   if i == 1
       hold on;
   end
end
hold off;
Legend = cell(length(SNR), 1);
for i = 1:length(SNR)
    Legend{i} = strcat('Signal to Noise Ratio: ', num2str(SNR(i)));
end
legend(Legend);

f6 = figure(6);
title('Unwrapped angle of signals');
ylabel('Angle [rad]');
xlabel('SNR [dB]');
grid on;
for i = 1:length(SNR)
    sig = A * exp(1j * (est(i, 1) * n' * T + mod(est(i, 2), pi)));
    plot(n, unwrap(angle(x(i, :))));
    if i == 1
        hold on;
    end
end
hold off;

Legend = cell(length(SNR), 1);
for i = 1:length(SNR)
    Legend{i} = strcat('Signal: ', num2str(i));
end
legend(Legend);
saveas(f6, 'figures/unwrapped_angles', 'epsc');

% H = [T * n', ones(length(n), 1)];
% C = A^2 / (2 * db2mag(SNR(length(SNR)))) * eye(N);
% est = (inv(H' * inv(C) * H) * H' * inv(C)) * unwrap(angle(x(length(SNR), :)))';
% est = BLUE(x(length(SNR), :), H, C);
% 
% w_0_est = est(1);
% phi_est = est(2);
% 
% x_gen = unwrap(angle(x(length(SNR), :)));
% x_est = unwrap(angle(A * exp(1j * (w_0_est * n' * T + phi_est))));
% 
% plot(n, x_gen, n, x_est);
% legend("x_gen", "x_est");

% plot(n, unwrap(angle(x(1, :))), n, unwrap(angle(x(2, :))), n, unwrap(angle(x(3, :))), n, unwrap(angle(x(4, :))), n, unwrap(angle(x(4, :))), n, unwrap(angle(x(6, :))));
% legend("SNR -10", "SNR 0" , "SNR 10" , "SNR 20" , "SNR 30", "SNR 40");


diff = zeros(length(SNR), N);
for snr = 1:length(SNR)
    for i = 1:N-1
        diff(snr, i) = angle(x(snr, i + 1)) - angle(x(snr, i));
    end
end

H_coloured = T * n';




