clear all
close all
clc



%% MAIN FUNCTION
% Notes Basics
k_list =  0:11;
abbs = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"];
names = ["ÕÙÔ", "ÕÙÔ#", "—Â", "—Â#", "ÃÈ", "÷·", "÷·#", "”ÔÎ", "”ÔÎ#", "À·", "À·#", "”È"];
oktava_list = 1:8;
f = zeros(length(oktava_list), length(k_list));
for ii = 1 : length(oktava_list)
    for jj = 1 : length(k_list)
        oktava = oktava_list(ii);
        k = k_list(jj);
        f(ii, jj) = nota_Hz(k, oktava);
    end
end

% Music Basics
Fs = 16000;
step = 1/Fs;
olokliro = 4;                       % Hz
miso = olokliro / 2;                % Hz
tetarto = miso / 2;                 % Hz
ogdoo = tetarto / 2;                % Hz
dekatoekto = ogdoo / 2;             % Hz

% Notes of the song and oktava selection
oktava = 2;
notes_names = ["”ÔÎ", "ÃÈ", "ÃÈ", "÷·", "—Â", "—Â", ...
               "ÕÙÔ", "—Â", "ÃÈ", "÷·", "”ÔÎ", "”ÔÎ", "”ÔÎ", ...
               "”ÔÎ", "ÃÈ", "ÃÈ", "÷·", "—Â", "—Â", ...
               "ÕÙÔ", "ÃÈ", "”ÔÎ", "”ÔÎ", "ÕÙÔ"];
           
LEN = length(notes_names);
notes_f = zeros(1, LEN);
notes_durations = ogdoo * ones(1, LEN);
stops_durations = dekatoekto/5 * ones(1, LEN-1);

% Redefine the stops
stops_durations(3) = ogdoo;
stops_durations(6) = ogdoo;
stops_durations(13) = ogdoo;
stops_durations(16) = ogdoo;
stops_durations(19) = ogdoo;

for n = 1 : LEN
    note_name = notes_names(n);
    k = find(names == note_name);
    notes_f(n) = f(oktava, k);
end




%% Play the songs
% ORIGINAL SONG
fprintf('<strong>Original</strong>\n\n');
offset_Hz_original = 0;
mult_original = 1;
original = create_song(notes_names, notes_f, notes_durations, ...
                   stops_durations, LEN, step, Fs, offset_Hz_original, mult_original, names, abbs);
sound(original, Fs);
pause(1.05 * sum([notes_durations stops_durations]));

% PARAFONO SONG
fprintf('\n<strong>Parafono</strong>\n');
offset_Hz_parafono = 50;
mult_parafono = 1;
disp("Mult = " + num2str(mult_parafono) + ", offset = " + num2str(offset_Hz_parafono) + " Hz");
display(' ');
display(' ');
parafono = create_song(notes_names, notes_f, notes_durations, ...
                   stops_durations, LEN, step, Fs, offset_Hz_parafono, mult_parafono, names, abbs);
sound(parafono, Fs);




%% Auxiliary Functions
function f_Hz = nota_Hz(k, oktava)
    f_Hz = 110 * 2^((k-9)/12) * 2^(oktava-1);
end

function song = create_song(notes_names, notes_f, notes_durations, ...
                            stops_durations, LEN, step, Fs, offset_Hz, mult, names, abbs)

    song = [];
    for n = 1 : LEN-1
        
        % Play the note
        t = 0 : step : notes_durations(n);
        f_Hz = notes_f(n);
        play = sin(2*pi* (mult*f_Hz + offset_Hz) *t);
        % Stop - pause, delay
        t = 0 : step : stops_durations(n);
        f_stop = 0;
        stop = sin(2*pi*f_stop*t);
        % Add to the song
        song = [song play stop];
        
        % Display the note
        note_name = notes_names(n);
        index = find(names == note_name);
        fprintf('%s <----> %.2f Hz (%.2f Hz deviation)\n', abbs(index), mult*f_Hz + offset_Hz, mult*f_Hz + offset_Hz - f_Hz);
        
    end

    % Add the last note
    t = 0 : step : notes_durations(n);
    f_Hz = notes_f(LEN);
    play = sin(2*pi* (mult*f_Hz + offset_Hz) *t);
    song = [song play];
    
    note_name = notes_names(LEN);
    index = find(names == note_name);
    fprintf('%s <----> %.2f Hz (%.2f Hz deviation)\n', abbs(index), mult*f_Hz + offset_Hz, mult*f_Hz + offset_Hz - f_Hz);
    
    
end