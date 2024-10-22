clc
clear all

% Data standard
Fs = 4000;
T = 1 / Fs;
LEN = 8000;
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
x = s + 0.75 * randn(size(t));

% 1. Define the frequency domain based on 'Fs'
f = Fs/LEN * (0:LEN/2);
% 2. Calculate the fft's
S = fft(s);
X = fft(x);
% 3. Two sided spectrums
% two_sided = P2 and single_sided = P1
P2_S = abs(S / LEN);
P1_S = P2_S(1 : 1+LEN/2);
P1_S(2:end-1) = 2 * P1_S(2:end-1);
P2_X = abs(X / LEN);
P1_X = P2_X(1 : 1+LEN/2);
P1_X(2:end-1) = 2 * P1_X(2:end-1);

% 4. Plots (time, two sided spectrum, single sided)
stop = LEN / 50;
figure();
title_all = "A_1=" + num2str(A1) + ", A_2=" + num2str(A2) + ...
          ", f_1=" + num2str(f1) + ", f_2=" + num2str(f2);
title(title_all);

subplot(3, 2, 1);
plot(1000*t(1:stop), s(1:stop));
title("Original in time (first samples)");
xlabel('t (milliseconds)')
ylabel('s(t)')
subplot(3, 2, 3);
plot(real(P2_S));
title("Two sided spectrum - original");
subplot(3, 2, 5);
plot(f, real(P1_S));
title("Single sided spectrum - original");
xlabel('f (Hz)')
ylabel('|P1_S(f)|');

subplot(3, 2, 2);
plot(1000*t(1:stop), x(1:stop));
title("Corrupted in time (first samples)");
xlabel('t (milliseconds)')
ylabel('x(t)')
subplot(3, 2, 4);
plot(real(P2_X));
title("Two sided spectrum - corrupted");
subplot(3, 2, 6);
plot(f, real(P1_X));
title("Single sided spectrum - corrupted");
xlabel('f (Hz)')
ylabel('|P1_X(f)|');

% Sound
sound(s, Fs);
delay(2.2);
sound(x, Fs);



% AUXILIARY FUNCTIONS
function y = delay(seconds)
    c = 0;
    for i = 1 : seconds * 10^9;
        c = c + 1;
    end
    y = c;
end

