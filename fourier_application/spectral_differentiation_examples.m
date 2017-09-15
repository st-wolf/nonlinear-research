%%
clc
clear

N = 24; h = 2*pi/N; x = h*(1:N)';
v = cos(x); v_hat = fft(v);

w_hat = 1i * [0:N/2-1 0 -N/2+1:-1]' .* v_hat;
w = real(ifft(w_hat));

plot(x, v, x, w);

figure; hold on
plot(real(w_hat));
plot(real(v_hat), 'Color', 'red');
 
%%
clc
clear

L = 2; N = 501; h = L/N; x = h*(1:N);
v = exp(cos(pi * x)); v_hat = fft(v);

% v = (x - 1) .^ 2 - 1; v_hat = fft(v);
% v = exp(sin(pi * x)); v_hat = fft(v);

w_hat = 1i * (2 * pi / L) * [0:(N-1)/2, -(N-1)/2:-1] .* v_hat;
w = real(ifft(w_hat));

plot(x, v, x, w);

figure; hold on
plot(real(w_hat));
plot(real(v_hat), 'Color', 'red');