clear all
clc

% Create random n-D spots
% "x" will contain 500 spots with "n" dimensions-coordinates
% Matrix will be Mxn, where M=500 spots with n coordinates

M = 500;
N_max = 50;
n_List = 1 : N_max;
MAX_list = zeros(1, N_max);
MIN_list = zeros(1, N_max);
LOG_list = zeros(1, N_max);

for k = 1 : length(n_List)
   
    n = n_List(k);
    spots = createSpots(M, n);
    distances = createDistances(M, n, spots);
    % distances(5, 7) is the distance between the 5-th and the 7-th
    % spots/rows in the spots matrix
    [MAX, MIN] = findMaxMin(M, n, distances);
    LOG = log10((MAX - MIN) / MIN);
    prettyPrint(M, n, MAX, MIN, LOG);
    MAX_list(k) = MAX;
    MIN_list(k) = MIN;
    LOG_list(k) = LOG;
    
end

plot(n_List, LOG_list);



% AUXILIARY FUNCTIONS
function spots = createSpots(M, n)
    spots = zeros(M, n);
    for i = 1 : M
        for j = 1 : n
            spots(i, j) = randn();
        end
    end
end

function distances = createDistances(M, n, spots)
    for elem1 = 1 : M
        for elem2 = 1 : M
            SUM = 0;
            for dim = 1 : n
                SUM = SUM + (spots(elem1, dim) - spots(elem2, dim))^2;
            end
            distances(elem1, elem2) = sqrt(SUM);
        end
    end
end

function [MAX, MIN] = findMaxMin(M, n, distances)
    MAX = max(max(distances));
    MIN = Inf;
    for i = 1 : M
        for j = 1 : M
            value = distances(i, j);
            if value > 0 && value < MIN
                MIN = value;
            end
        end
    end
end

function prettyPrint(M, n, MAX, MIN, LOG)
    display('*****************************************');
    disp("M = " + num2str(M) + " " + num2str(n)+ "-D spots");
    disp("Max = " + num2str(MAX)+ ", MIN = " + num2str(MIN) + ", LOG = " + num2str(LOG));
    display('*****************************************');
    display(' ');
end