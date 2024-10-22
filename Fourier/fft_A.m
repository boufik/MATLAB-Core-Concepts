clear all
close all
clc

%% A1. Basics
Fs = 500;                   % Sampling Frequency   [Hz]
duration = 2;               % Signal time duration [sec]
Ts = 1 / Fs;                % Sampling period      [sec]

t = 0 : Ts : duration - Ts;
N = length(t);              % Num of samples 
    

%% A2. Signal parameters
A1 = 3;     f1 = 30;    phase1 = 0.6;
A2 = 2;     f2 = 45;    phase2 = -0.8;
A3 = 1;     f3 = 70;    phase3 = 2;

s1 = A1*cos(2*pi*f1*t + phase1);
s2 = A2*cos(2*pi*f2*t + phase2);
s3 = A3*cos(2*pi*f3*t + phase3);
s = s1 + s2 + s3;


%% A3. Plotting in time-domain
figure();
plot(t, s);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Time-domain plot');


%% A4. Fast Fourier Transform - Samples
S = fft(s);                     % Also 'N' samples
figure();
plot(abs(S));
xlabel('Samples');
ylabel('Magnitude');
title('Samples-domain plot');


%% A5. Convert samples into frequencies by keeping the LHS of plot
S_oneside = S(1 : N/2);
f = Fs * (0 : N/2 - 1) / N;
S_meg = abs(S_oneside) / (N/2);

figure();
plot(f, S_meg);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Frequency-domain plot (Sharp Spikes)');


%% A6. Discover the sinus phases
ph1 = angle(S_oneside(f1*duration+1))
ph2 = angle(S_oneside(f2*duration+1))
ph3 = angle(S_oneside(f3*duration+1))


