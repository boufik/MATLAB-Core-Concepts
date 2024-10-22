clc
clear all

% Time
[d, Fs] = audioread('kanonikosima.mp3');
subplot(2, 1, 1);
plot(d);
xlabel('Time');
title('Time domain');

% Frequency
D = fft(d(:, 1));
% 3. Two sided spectrums, two_sided = P2 and single_sided = P1
LEN = size(d, 1);
P2_S = abs(D / LEN);
P1_S = P2_S(1 : 1+LEN/2);
P1_S(2:end-1) = 2 * P1_S(2:end-1);
P2_X = abs(D / LEN);
P1_X = P2_X(1 : 1+LEN/2);
P1_X(2:end-1) = 2 * P1_X(2:end-1);
subplot(2, 1, 2);
plot(P1_X);
xlabel('Frequency (Hz)');
title('One sided spectrum');