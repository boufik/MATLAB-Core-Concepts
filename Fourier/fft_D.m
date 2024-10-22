clear all
close all
clc

%% D1. New duration != integer number of signal cycles
Fs = 500;                   % Sampling Frequency   [Hz]
duration = 0.64;            % Signal time duration [sec]
Ts = 1 / Fs;                % Sampling period      [sec]
t = 0 : Ts : duration - Ts;
N = length(t);              % Num of samples 
    

%% D2. Signal parameters - Windowing BEFORE padding
A1 = 3;     f1 = 30;    phase1 = 0.6;
A2 = 2;     f2 = 45;    phase2 = -0.8;
A3 = 1;     f3 = 70;    phase3 = 2;
s1 = A1*cos(2*pi*f1*t + phase1);
s2 = A2*cos(2*pi*f2*t + phase2);
s3 = A3*cos(2*pi*f3*t + phase3);
s = s1 + s2 + s3;

s = s .* hanning(N)';
figure();
plot(t, s);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Time-domain plot - Signal after windowing');




%% C3. Plotting in time-domain - Zero-padding signal 's'
s = [s zeros(1, 2000)];
t_extra = duration : Ts : duration + 1999*Ts;
t = [t t_extra];

figure();
plot(t, s);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Time-domain plot - Windowing and padding');


%% C4. Fast Fourier Transform - Samples
S = fft(s);                     % Also 'N' samples
figure();
plot(abs(S));
xlabel('Samples');
ylabel('Magnitude');
title('Samples-domain plot');


%% C5. Convert samples into frequencies by keeping the LHS of plot
% Side lobes will disappear after having windowed and padded our signal
N2 = length(s);
S_oneside = S(1 : N2/2);
f = Fs * (0 : N2/2 - 1) / N2;
S_meg = abs(S_oneside) / (N2/2);
figure();
plot(f, S_meg);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Frequency-domain plot (Side lobes disappeared with windowing and padding)');