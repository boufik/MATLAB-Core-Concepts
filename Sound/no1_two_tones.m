clc
clear all

% Data standard
Fs = 4410;
T = 1 / Fs;
LEN = 1000;
t = (0 : LEN-1) * T;
% Selecting data
A1 = 0.7;
A2 = 0.9;
f1 = 294;
f2 = 440;
s1 = A1 * cos(2*pi*f1*t);
s2 = A2 * cos(2*pi*f2*t);
s = s1 + s2;
S1 = fft(s1);
S2 = fft(s2);
S = fft(s);

% 1. Figure 1 - Time and frequency plots along with sound
figure(1);
title("Time and frequency domain");

sound(s1, Fs);
delay(1);
title1 = num2str(f1) + " Hz";
subplot(3, 2, 1);
plot(s1);
title(title1 + " in time");
subplot(3, 2, 2);
plot(real(S1));
title(title1 + " in frequency");

sound(s2, Fs);
delay(1);
title2 = num2str(f2) + " Hz";
subplot(3, 2, 3);
plot(s2);
title(title2 + " in time");
subplot(3, 2, 4);
plot(real(S2));
title(title2 + " in frequency");

sound(s, Fs);
delay(1);
title3 = "Two tones together";
subplot(3, 2, 5);
plot(s);
title(title3 + " in time");
subplot(3, 2, 6);
plot(real(S));
title(title3 + " in frequency");





% AUXILIARY FUNCTIONS
function y = delay(seconds)
    c = 0;
    for i = 1 : seconds * 10^9;
        c = c + 1;
    end
    y = c;
end

