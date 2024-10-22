clear all
clc
 
 
% Main Function
% Standard data
M = 24;
N = 13;
Pg = 1;
SNR = 10;
thetaList = [];
counter = 0;
for angle = 30:10:150
    counter = counter + 1;
    thetaList(counter) = angle * pi/180;
end
THETA_LIST = thetaList * 180/pi
run(M, N, Pg, SNR, thetaList);
 
 
% ===============================================================
% ===============================================================
% ======================== Functions ============================
% ===============================================================
% ===============================================================
 
% Function 1 - Steering vector
function a = steeringVector(M, theta)
    a = [];    
    for m = 1 : M
        a(m) = exp(i * pi * (m-1) * cos(theta));
    end
    a = a';
    a;
end
 
% Function 2 - Steering Matrix
function A = steeringMatrix(M, thetaList)
    A = [];
    for index = 1 : length(thetaList)
        A(:, index) = steeringVector(M, thetaList(index));
    end
    A;
end
 
% Function 3 - DoA: LP = Linear Prediction
function spectrum = draw_LP_Spectrum(M, N, A,Pg, P_noise, thetaList)
    Rgg = Pg * eye(N);
    Rnn = P_noise * eye(M);
    Rxx = A * Rgg * A' + Rnn;
    % Draw χωρικό φάσμα ισχύος του εκτιμητή
    % We suppose 1st element as the reference
    e1 = zeros(1, M);
    e1(1) = 1;
    e1 = e1';            % Now, its a column-vector with size = 24x1
    nominator = e1' * inv(Rxx) * e1;
    % For denominator, I need to make a vector "ad" which refers to a random
    % angle signal
    counter = 0;
    angles = [];            % x axis
    spectrum = [];          % y axis
    for angle = 0:0.1:180
        counter = counter + 1;
        angles(counter) = angle;        % In degrees
        angleRad = angle * pi/180;      % In rad
        ad = steeringMatrix(M, angleRad);
        % For calculations, I will need rad angles
        denominator = (norm(e1' * inv(Rxx) * ad))^2;
        spectrum(counter) = nominator / denominator;
    end
    plot(angles, spectrum);
    myTitle = "";
    THETA_LIST = thetaList * 180/pi;
    for i = 1:length(THETA_LIST)
        t = "θ" + int2str(i) + "=" + num2str(THETA_LIST(i));
        if i < length(THETA_LIST)
            t = t + ", ";
        end
        myTitle = myTitle + t;
    end
    title(myTitle);
    xlabel("θ in degrees");
    ylabel("10log(P / Pmax)");
    spectrum;
end
 
% Function 4 - Run with data
function run(M, N, Pg, SNR, thetaList)
    SNR_noDB = 10^(SNR/10);
    P_noise = Pg / SNR_noDB;
    A = steeringMatrix(M, thetaList);        % Size = 24x13
    spectrum = draw_LP_Spectrum(M, N, A, Pg, P_noise, thetaList);
end
