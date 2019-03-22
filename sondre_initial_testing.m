clc; clear;


F_s = 10^6;
T = 1/F_s;
f_0 = 10^5;
w_0 = 2*pi*f_0;
phi = pi/8;
A = 1;

%N samples from n = n_0 through n = n_0 + N-1
%in b it says N = 513, n_0 = -256
N = 513;
P = (N*(N-1)/2);
Q = (N*(N-1)*(2*N-1))/6;

n_0 = -P/N; n_N = n_0+N;

v = (rand([N 1]) + 1i*rand([N 1])) .*2 -1;

x = A*exp(1i*(w_0*(n_0:n_N)'*T + phi + v));

unwrapped = unwrap(angle(x));
%% Plotting

figure(1)
plot(1:N,angle(x), 1:N, unwrapped)
legend('angle(x)', 'unwrapped')
%stem(x)
title('x[n]')





