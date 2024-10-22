clear all
close all
clc

%% 1a. Create a periodic signal as the sum of 2 sinus waves
syms t n
STOP = 24;
tt = linspace(0, STOP, 1201);
T1 = 3;     T2 = 2;     A1 = 2.5;   A2= 1.5;
yy = A1*sin(2*pi*tt / T1) + A2*sin(2*pi*tt / T2);

figure();
plot(tt, yy);
title("Fourier series");
xlabel("t");
hold on


%% 1b. Find the Fourier series coefficients
T = lcm(T1, T2);  f0 = 1/T;   w0  = 2*pi*f0;
N_LIMITS = [10 20 50];
titles = ["Original pulse"];

for index = 1 : length(N_LIMITS)
    
    display('************************************************');
    N_LIMIT = N_LIMITS(index);
    n = 1 : N_LIMIT;
    fprintf('N = %d\n', N_LIMIT);
    
    y = A1*sin(2*pi*t / T1) + A2*sin(2*pi*t / T2);
    a0 = (1/T) * int(y, t, 0, 2)
    an = (2/T) * int(y*cos(n*w0*t), t, 0, 2)
    bn = (2/T) * int(y*sin(n*w0*t), t, 0, 2)

    y_appr = @(t) a0;
    for i = 1 : N_LIMIT
        an_ = an(i);
        bn_ = bn(i);
        n_ = n(i);
        y_appr = y_appr + an_ * cos(n_*w0*t) + bn_ * sin(n_*w0*t);
    end
    fplot(t, y_appr, [0 STOP]);
    hold on
    titles(end+1) = num2str(N_LIMIT) + "-terms Fourier series";
    fprintf('************************************************\n\n');
end

legend(titles);


