clear all
clc
 
 
% Main Function
% Standard data
M = 24;
Pg = 1;    
minAngle = 30*pi/180;
maxAngle = 150*pi/180;
stepDegrees = 0.1;
% Delta and SNR
%{
deltaList = [6, 8, 10, 12, 14, 16] * pi/180;
SNRlist = [0, 5, 10, 15, 20, 25, 30];
for i = 1:length(deltaList)
    delta = deltaList(i)
    for j = 1:length(SNRlist)
        SNR = SNRlist(j)
        run(M, delta, SNR, Pg, minAngle, maxAngle);
    end
end
%}
SNR = 20;
delta = 10*pi/180;
run(M, delta, SNR, Pg, minAngle, maxAngle, stepDegrees);
 
 
 
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
 
% Function 3 - Find the maximum angle0 (minimum is 30 degrees always)
function maxAngle0 = findMaxAngle0(maxAngle, delta)
    maxAngle0 = maxAngle - 5 * delta;
end
 
% Function 4 - Create arxiki exada
function arxikiExada = createArxikixada(angle0, delta)
    arxikiExada = [angle0, angle0 + delta, angle0 + 2*delta, angle0 + 3*delta, angle0 + 4*delta, angle0 + 5*delta];
end
 
% Function 5 - Having a single exada, I will create the new 6 exades
% icluding the 6 differents combinations for theta0
function exades = createExades(arxikiExada)
    exada1 = arxikiExada;
    exada2 = [arxikiExada(2) arxikiExada([1, 3:end])];
    exada3 = [arxikiExada(3) arxikiExada([1:2, 4:end])];
    exada4 = [arxikiExada(4) arxikiExada([1:3, 5:end])];
    exada5 = [arxikiExada(5) arxikiExada([1:4, end])];
    exada6 = [arxikiExada(6) arxikiExada([1:5])];
    % Exada_i = a row vector with dimensions 1x6
    % My output will be an matrix 6x6, each row is a new exada
    exades = [exada1; exada2; exada3; exada4; exada5; exada6];
end
 
 
% Function 6 - Calculate the NSB weights given an exada = thetaList 
function weights = NSBweights(M, P_noise, exada)
    A = steeringMatrix(M, exada);
    tempMatrix = A' * A + P_noise * eye(6);
    invTemp = inv(tempMatrix);
    e1 = [1; 0; 0; 0; 0; 0];
    weights = A * invTemp * e1
end
 
 
% Function 7 - Radiation diagram and calculation of SINR
function AF = af(M, weights, exada, stepDegrees)
    counter = 0;
    angles = [];
    AF = [];
    for theta = 0: stepDegrees*pi/180: pi
        counter = counter + 1;
        angles(counter) = theta * 180/pi;
        AF_complex = weights' * steeringVector(M, theta);
        AF(counter) = abs(AF_complex);
    end
    % Info in the title
    theta0 = exada(1)*180/pi;
    theta1 = exada(2)*180/pi;
    theta2 = exada(3)*180/pi;
    theta3 = exada(4)*180/pi;
    theta4 = exada(5)*180/pi;
    theta5 = exada(6)*180/pi;
    titleInfo = "è0 = " + int2str(theta0) + ", è1 = " +...
        int2str(theta1) +  ", è2 = " + int2str(theta2) +...
        ", è3 = " + int2str(theta3) + ", è4 = " +...
        int2str(theta4) + ", è5 = " + int2str(theta5); 
    plot(angles, AF);
    title(titleInfo);
    xlabel("è (in degrees)");
    ylabel("|AF(è)|");
    AF;
end
 
% Function 8 - Calculate SINR
function SINR = calculateSINR(M, weights, exada, Pg, P_noise)
    % I need to find Rdd (theta0) and Ruu (theta1 to theta5)
    % Rdd = Psd * ad * ad'
    desiredTheta = exada(1);
    ad = steeringVector(M, desiredTheta);
    undesiredAngles = exada([2:6]);
    Ai = steeringMatrix(M, undesiredAngles);
    
    Rdd = Pg * ad * ad';
    Rgigi = Pg * eye(5);        % Because, interferences have same power
    Ruu = Ai * Rgigi * Ai' + P_noise * eye(M);
    
    Pyd = weights' * Rdd * weights;
    Pyu = weights' * Ruu * weights;
    SINR_noDB = abs(floor(Pyd / Pyu));
    SINR = 10 * log10(SINR_noDB);
end
 
% Function 9 - Calculate maximum side lobes level
function SLL_and_theta0_divergence = calculateSLL(AF, exada, stepDegrees)
    % TO find SLL, first I will find all the topic max
    topicMaxList = [];
    counterMax = 0;
    for i = 2:length(AF)-1
        if AF(i) > AF(i-1)&& AF(i) > AF(i+1)
            counterMax = counterMax+1;
            topicMaxList(counterMax) = AF(i);
        end
    end
    % After for-loop, I will scan topicMaxList to find the 2 max values
    % 1st one is the height of main lobe and the 2nd one is for the
    % most powerful side lobe
    mainLobe = max(topicMaxList)
    index1 = find(topicMaxList == mainLobe, 1);
    % Delete element in AF
    topicNew = topicMaxList([1:index1-1, index1+1:end]);
    sideLobe = max(topicNew)
    index2 = find(topicNew == sideLobe, 1);
    SLL = 10 * log10(mainLobe / sideLobe);
    % After that I have to find the theta in which I have the mainLobe
    % value scanning the AF array
    index = find(AF == mainLobe, 1);
    thetaDegrees = (index - 1) * stepDegrees;   % This is not the desired one
    desired = exada(1) * 180/pi;                % This is the desired one
    theta0_divergence = abs(desired - thetaDegrees);
    SLL_and_theta0_divergence = [SLL theta0_divergence];  % Return 2 values
    % Return mainLobe height and position
end
 
 
% Function 10 - Calculate divergences
function divergences = calculateDivergences(AF, exada, stepDegrees)
    undesired = exada([2:6]) * 180/pi      % Degrees (1x5 array)
    epsilon = 0.001;                % Upper bound (lower = 0)
    counterMin = 0;
    zeroThetaList = [];
    for i = 1:length(AF)
        if AF(i) <= epsilon
            counterMin = counterMin + 1;
            zeroThetaList(counterMin) = i * stepDegrees;
            % index i indicates degree
        end
    end
    zeroThetaList;
    % Now, my list contains the angles with the lowest possible values of AF
    % I have to find the distance between them and my undesired thetas
    divergences = [];
    for i = 1:length(undesired)
        undesiredTheta = undesired(i);
        minDistance = 1000;
        for j = 1:length(zeroThetaList)
            zeroTheta = zeroThetaList(j);
            if abs(undesiredTheta - zeroTheta) <= minDistance
                minDistance = abs(undesiredTheta - zeroTheta);
                % It is about this specific undesiredTheta
            end
        end
        divergences(i) = minDistance;
    end
    divergences
end
 
 
% Function 11 - Statistics into a given vector
function statistics = calculateStatistics(myList)
    maxElement = max(myList);
    minElement = min(myList);
    SUM = sum(myList);
    mean = SUM / length(myList);
    % STD
    SUM2 = 0;
    for i = 1:length(myList)
        SUM2 = SUM2 + (myList(i) - mean)^2;
    end
    std = sqrt(SUM2 / length(myList));
    statistics = [minElement maxElement mean std];
end
 
 
 
 
% ===============================================================
% ===============================================================
% Function 12 - Run with data
function run(M, delta, SNR, Pg, minAngle, maxAngle, stepDegrees)
    % Example for angle delta = 10 moires
    % counter = 1  <----> 30 and 40,50,60,70,80
    % counter = 7  <----> 40 and 50,60,70,80,90
    % counter = 13 <----> 50 and 60,70,80,90,100
    % counter = 19 <----> 60 and 70,80,90,100,110
    counter = 0;
    SNR_noDB = 10 ^ (SNR / 10);
    P_noise = Pg / SNR_noDB;
    % Testing 3
    minAngle0 = minAngle;
    maxAngle0 = findMaxAngle0(maxAngle, delta);
    % For statistic purposes
    theta0_divergencesList = [];      % pe 71 elements
    divergencesList = [];             % pe 355 = 71 * 5 elements
    SINRlist = [];                    % 71 elements
    SLLlist = [];                     % 71 elements
    % Testing 4, 5, 6, 7, 8, 9
    % Variable "i" will describe EACH UNIQUE SET OF 6 FIRST ANGLES!!!!
    for angle0 = minAngle0 : pi/180 : maxAngle0   % We work with rads, so 1 degree = pi/180 rad
        display(angle0 * 180/pi);   
        arxikiExada = createArxikixada(angle0, delta);    % Connected with "i"
        exades = createExades(arxikiExada);
        for row = 1:6
            counter = counter + 1
            exada = exades(row, :);
            display(exada * 180/pi);
            
            weights = NSBweights(M, P_noise, exada);
            AF = af(M, weights, exada, stepDegrees);
            SINR = calculateSINR(M, weights, exada, Pg, P_noise);
            SLL_and_theta0_divergence = calculateSLL(AF, exada, stepDegrees);
            SLL = SLL_and_theta0_divergence(1)
            theta0_divergence = SLL_and_theta0_divergence(2)
            divergences = calculateDivergences(AF, exada, stepDegrees);
   
            % For statistic purposes
            theta0_divergencesList(counter) = theta0_divergence;
            D = length(divergences);
            LENGTH = length(divergencesList);
            for j = 1:D
                divergencesList(LENGTH + j) = divergences(j);
            end
            SINRlist(counter) = SINR;
            SLLlist(counter) = SLL;
            
            display('******************************************************');
            display(' ');
        end
    end
    theta0_divergencesList;
    divergencesList;
    SINRlist;
    SLLlist;
    display(' ');
    display(' ');
    display('******************************************************');
    display(strcat('Delta=', int2str(delta*180/pi), ' and SNR=', int2str(SNR)));
    display(' ');
    theta0_divergences_stats = calculateStatistics(theta0_divergencesList)
    divergences_stats = calculateStatistics(divergencesList)
    SINR_stats = calculateStatistics(SINRlist)
    SLL_stats = calculateStatistics(SLLlist)
end
% ===============================================================
% ===============================================================
 
 
% Gia ta prin, esvisa
% delta (3), maxExades (3), arxikiExada (4), exades (5), AF(7)