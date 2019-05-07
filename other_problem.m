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
v = zeros(length(SNR), N);
v = v;
diff = zeros(length(SNR), N - 1);
var = (A^2 / 2) ./ db2mag(SNR);
est = zeros(length(SNR), 2);
H = T * ones(N-1, 1);
D_base = diag([ones(1, N-1), 0]);
D_base = D_base(1:N-1, 1:N);
D = circshift(D_base, 1, 2) - D_base;

%hpFilter = designfilt('bandpassiir', 'FilterOrder', 10, ...
%    'HalfPowerFrequency1', 80, 'HalfPowerFrequency2', 120, 'SampleRate', 500);

for i = 1:length(SNR)
    x(i, :) = gen_signal(w_0, n, A, T, phi, 0, sqrt(var(i)));
    v(i, :) = gen_noise(N, 0, sqrt(var(i)));
    %v(i,:) = filter(hpFilter, v(i,:)); %coloring noise
    y = angle(x(i, :))';
    for j = 1:N-1
        %diff(i, j) = y(j + 1) - y(j);
        diff(i, j) = w_0*T + v(j+1) - v(j);
    end
    
    C = D * (var(i)) * D';
    % est(i, 1) = BLUE_c(diff(i, :)', H, C);
    est(i, 1) = abs(BLUE_c(diff(i, :)', H, C));
    
    %fourier = F(N, est(i, 1), T);
    res = zeros(N, 1);
    for n_ = 1:N
        res(n_) = x(i, n_) * exp(-1j * est(i, 1) * (n_ - 1)*T);
    end
    
    fourier = mean(res);
    %fourier = F(x(i, :), N, w_0, T);
    est(i, 2) = angle(exp(-1j*est(i,1)*n_0*T)*fourier); %something wierd goin on here
        %It is the line stated in the assigment text, but it's not giving
        %the correct result.
    
    % C_inv = inv(var(i) * (D * D'));
    % est(i, 1) = inv(H' * C_inv * H) * (H' * C_inv * (D * angle(x(i, :)')));
    %inv(H' * inv(var(i) * D * D') * H) * H' * inv(var(i) * D * D') * D * angle(x(i, :)');
    %BLUE_c(D * angle(x(i, :)'), H, D * var(i) * eye(N) * D');
    %est(i, 2) = (1 / N) * sum(angle(x(i, :) - est(i, 1) * n));
end

CRLB_omega = (12 / (A^2 * T^2 * N * (N^2 - 1))) .* var;

%Filter here



wfft = fft(v(6,:));     
f = (0:length(wfft)-1)*50/length(wfft);

plot(f,abs(wfft))

