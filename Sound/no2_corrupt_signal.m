clc
clear all

% Data standard
Fs = 5000;
T = 1 / Fs;
LEN = 10000;
t = (0 : LEN-1) * T;
% Selecting data
A1 = 0.7;
A2 = 0.9;
f1 = 294;
f2 = 440;
s1 = A1 * cos(2*pi*f1*t);
s2 = A2 * cos(2*pi*f2*t);
s = s1 + s2;
% Corrupt my signal s ---> create x
x = s + 0.35 * randn(size(t));

% Display the first time samples
S = fft(s);
X = fft(x);
stop = LEN / 50;         % Show only 2% of all the samples
figure();

subplot(2, 2, 1);
plot(1000*t(1:stop), s(1:stop));
title("Original in time (first samples)");
xlabel('t (milliseconds)')
ylabel('s(t)')
subplot(2, 2, 2);
plot(real(S));
title("Original in frequencies (all samples)");

subplot(2, 2, 3);
plot(1000*t(1:stop), x(1:stop));
title("Corrupted in time (first samples)");
xlabel('t (milliseconds)')
ylabel('x(t)')
subplot(2, 2, 4);
plot(real(X));
title("Corrupted in frequencies (all samples)");

% Sound
sound(s, Fs);
display("Playing the ORIGINAL sound");
delay(2);
sound(x, Fs)
display("Playing the CORRUPTED sound");



% AUXILIARY FUNCTIONS
function y = delay(seconds)
    c = 0;
    for i = 1 : seconds * 10^9;
        c = c + 1;
    end
    y = c;
end

