clc
clear all




% MAIN FUNCTION
figure();
filenames = ["aa", "ee", "oo", "ii"];
for i = 1:length(filenames)
    filename = filenames(i);
    y = plots_sound(filename + ".mp3", i);
    delay(1);
    
end

% y = [1, 10, 2, 1, 3, 4, 2, 10, 11, 1, 3, 5, 6, 5];



% Function 1
function P1 = plots_sound(filename, i)
    [s, Fs] = audioread(filename);
    LEN = length(s(:, 1));
    t = (0:LEN-1) / Fs;

    % Create the single sided spectrum of time signal s(t)
    stop = floor(LEN / 50);
    S = fft(s);
    f = Fs/LEN * (0:LEN/2);
    P2 = abs(S / LEN);                  % Two sided P2
    P1 = P2(1 : 1+LEN/2);               % Single sided P1
    P1(2:end-1) = 2*P1(2:end-1);

    subplot(2, 2, i);
    plot(f, P1);
    title(filename);
    xlabel("f (Hz)");
    ylabel("|S(f)|");

    % Sound
    sound(s, Fs);
    topic_maxes = find_topic_maxes(P1);
    freq = topic_maxes(1, :);
    ampl = topic_maxes(2, :);
    [AMPL, INDEX] = max(ampl);
    FREQ = freq(INDEX);
    FREQ = FREQ * Fs / LEN;
    disp("Freq of " + filename + " ---> " + num2str(FREQ) + " Hz (" + num2str(AMPL) + ")");
    
end



% Function 2
function c = delay(n)
    c = 0;     
    for i = 1 : n * 10^9
        c = c + 1;
    end
end



% Function 3
function topic_maxes = find_topic_maxes(y)
    
    topic_maxes_pos = [];
    topic_maxes_val = [];
    for i = 2 : length(y) - 1
        if y(i) >= y(i-1) && y(i) >= y(i+1)
            topic_maxes_pos(end + 1) = i;
            topic_maxes_val(end + 1) = y(i);
        end
    end
    topic_maxes = [topic_maxes_pos ; topic_maxes_val];

end
