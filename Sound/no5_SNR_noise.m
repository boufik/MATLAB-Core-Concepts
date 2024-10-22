clc
clear all

[d, Fs] = audioread('kanonikosima.mp3');
duration = size(d, 1) / Fs;
sound(d, Fs);
subplot(3, 2, 1);
plot(d);
xlabel('Time');
ylabel('Sound');
title('Original sound');
display('Original sound');
delay(duration);
avg1 = mean(d(:, 1));
power1 = avg1 ^ 2;
SNR_dB_list = [-30, -40, -50, -60, -70];

for k = 1:length(SNR_dB_list)
    SNR_dB = SNR_dB_list(k);
    SNR = 10 ^ (SNR_dB / 10);
    power2 = power1 / SNR;
    avg2 = sqrt(power2);
    noise = zeros(size(d, 1), size(d, 2));
    for i = 1 : size(d, 1)
        for j = 1 : size(d, 2)
            noise(i, j) =  avg2 * randn();
        end
    end
    d_noise = d + noise;
    sound(d_noise, Fs);
    subplot(3, 2, k+1);
    plot(d_noise);
    xlabel('Time');
    ylabel('Sound');
    message = "Sound with SNR = " + num2str(SNR_dB) + " dB";
    title(message);
    disp(message);
    delay(duration);
end    



function y = delay(duration)
    c = 0;
    for i = 1 : duration * 1.15 * 10^9
        c = c + i;
    end
end