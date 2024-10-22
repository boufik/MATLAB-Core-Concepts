clear all
clc
 
 
% Main Function
% Standard data
M = 24;
N = 13;
Pg = 1;
SNR = 10;
stepDegrees = 0.1;
SNR_noDB = 10^(SNR / 10);
P_noise = Pg / SNR_noDB;
centerDegrees = 90;
counter = 0;
deltaDegreesList = [];
for i = 5:-0.01:4
    counter = counter + 1;
    deltaDegreesList(counter) = i;
end

% Scanning my list for all the values of delta (between 4 and 5 degrees)
for i = 1:length(deltaDegreesList)
    display('*********************************************');
    deltaDegrees = deltaDegreesList(i)
    % if deltaDegrees == 4.21
    delta = deltaDegrees * pi/180;
    thetaList = createThetaList(centerDegrees, deltaDegrees);
    THETA_LIST = thetaList * 180/pi
    A = steeringMatrix(M, thetaList);
    spectrum = draw_LP_Spectrum(M, N, A, THETA_LIST, Pg, P_noise, centerDegrees, deltaDegrees, stepDegrees);
    topicMax = findTopicMax(spectrum, stepDegrees);
    LENGTH = topicMax(1)
    topicMaxList = topicMax([2:LENGTH+1]);
    anglesList = topicMax([LENGTH+2:end])
    
    % Now, we have the angles of topic max, but sometimes their number is not
    % 13, so it's wrong
    if LENGTH ~= 13
        display(strcat('LP method is not efficient for delta=', num2str(deltaDegrees),' degrees'));
    end
    display(' ');
    % end
end

 
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
function spectrum = draw_LP_Spectrum(M, N, A, THETA_LIST, Pg, P_noise, centerDegrees, deltaDegrees, stepDegrees)
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
    for angle = 0:stepDegrees:180
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
    for i = 1:length(THETA_LIST)
        t = "θ" + int2str(i) + "=" + int2str(THETA_LIST(i));
        if i < length(THETA_LIST)
            t = t + ", ";
        end
        myTitle = myTitle + t;
    end
    title(myTitle);
    xlabel("θ in degrees, delta = " + num2str(deltaDegrees));
    ylabel("10log(P / Pmax)");
    spectrum;
end

% Function 4 - createThetaList(center, delta) - output in rad
function thetaList = createThetaList(centerDegrees, deltaDegrees)
    thetaList = [];
    counter = 0;
    for i = -6:6
        counter = counter + 1;
        angle = centerDegrees + i * deltaDegrees;     % Degrees
        theta = angle * pi/180;         % Rad
        thetaList(counter) = theta;
    end
    thetaList;
end


% Function 5 - Deleting zeros from a list
function newList = deleteZerosFrom(myList)
    counter = 0;    
    newList = [];
    for i = 1:length(myList)
        if myList(i) ~= 0
            counter = counter + 1;
            newList(counter) = myList(i);
        end
    end
    newList;
end

% Function 6 - Find the topic max values in my spectrum
function topicMax = findTopicMax(spectrum, stepDegrees)
    lowerBound = 20;
    counter = 0;
    topicMaxList = [];
    anglesList = [];
    % If my sample is below this value, I continue with the next sample
    for i = 2:length(spectrum)-1
        value = spectrum(i);
        if value >= lowerBound
            counter = counter + 1;
            if value > spectrum(i-1) && value > spectrum(i+1)
                topicMaxList(counter) = value;
                index = i;
                angle = (index - 1) * 0.1;
                anglesList(counter) = angle;
            end
        end
    end
    % Now, I have a list with all my topic max and the positions
    % I will return a big vector with 3 returned values in row
    % 1st_element = length = number of topic max values
    % 2nd_element = a list with all the topic max values (13 elememnts for example)
    % 3rd_element = a list with all the angles of topic max (13 elememnts for example)
    
    % Sometimes my 2 lists have some zeros (depends on lower bound)
    % I have to delete them
    topicMaxList = deleteZerosFrom(topicMaxList);
    anglesList = deleteZerosFrom(anglesList);
    LENGTH = length(topicMaxList);
    topicMax = [LENGTH topicMaxList anglesList];
end


