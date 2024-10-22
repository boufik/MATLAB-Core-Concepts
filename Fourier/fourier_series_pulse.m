clear all
close all
clc

%% 1a. Create a pulse with duty cycle = 0.5
syms t n
tt = linspace(0, 12, 1201);
yy = zeros(1, length(tt));
for i = 1 : length(tt)
    t_ = tt(i);
    if (0<=t_ && t_<=2) || (4<=t_ && t_<=6) || (8<=t_ && t_<=10)
        yy(i) = 1;
    end
end

figure();
plot(tt, yy);
title("Fourier series");
xlabel("t");
hold on


%% 1b. Find the Fourier series coefficients
T = 4;  f0 = 1/T;   w0  = 2*pi*f0;
N_LIMITS = [10 20 50];
titles = ["Original pulse"];

for index = 1 : length(N_LIMITS)
    
    display('************************************************');
    N_LIMIT = N_LIMITS(index);
    n = 1 : N_LIMIT;
    fprintf('N = %d\n', N_LIMIT);

    a0 = (1/T) * int(1, t, 0, 2)
    an = (2/T) * int(1*cos(n*w0*t), t, 0, 2)
    bn = (2/T) * int(1*sin(n*w0*t), t, 0, 2)

    y_appr = @(t) a0;
    for i = 1 : N_LIMIT
        an_ = an(i);
        bn_ = bn(i);
        n_ = n(i);
        y_appr = y_appr + an_ * cos(n_*w0*t) + bn_ * sin(n_*w0*t);
    end
    fplot(t, y_appr, [0 12]);
    hold on
    titles(end+1) = num2str(N_LIMIT) + "-terms Fourier series";
    fprintf('************************************************\n\n');
end

legend(titles);


