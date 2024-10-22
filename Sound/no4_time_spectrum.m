clc
clear all

% Declare time
Fs = 6000;                  % Max frequency of signal must be lower than Fs/2
LEN = 10000;
t = (0:LEN-1) / Fs;
% Ranges - Define the expression of signal s(t) - find Amplitudes and
% Frequencies
f_ranges = 200:1000;
A_ranges = 0.25:0.01:0.99;
NUM = 4;
f_list = f_ranges(randi(length(f_ranges), 1, NUM))
A_list = A_ranges(randi(length(A_ranges), 1, NUM))
s = 0;
for i = 1:NUM
    A = A_list(i);
    f = f_list(i);
    s = s + A * cos(2*pi*f*t);
end


% Create the single sided spectrum of time signal s(t)
stop = LEN / 50;
S = fft(s);
f = Fs/LEN * (0:LEN/2);
P2 = abs(S / LEN);                  % Two sided P2
P1 = P2(1 : 1+LEN/2);               % Single sided P1
P1(2:end-1) = 2*P1(2:end-1);

subplot(2, 1, 1);
plot(1000*t(1:stop), s(1:stop));
title("Time domain of 4 tones(first samples)");
xlabel("t (ms)");
ylabel("s(t) = 4 tones");

subplot(2, 1, 2);
plot(f, P1);
title("Single sided spectrum");
xlabel("f (Hz)");
ylabel("|P_1(f)|");

% Sound
sound(s, Fs);
