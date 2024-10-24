clear all
clc

% Data
N = 10^3;
x = 1 : N;
e = exp(1);

% Bagging Rounds + Original Data
MAX_ROUND = 5;
rounds_list = 0 : MAX_ROUND;
un_sim_list = zeros(1, MAX_ROUND+1);
un_exp_list = zeros(1, MAX_ROUND+1);

for i = 1 : length(rounds_list)
    round = rounds_list(i);
    if round == 0
        un_sim_list(round+1) = N;
        un_exp_list(round+1) = N;
    else
        x = BAGGING(x);
        un_sim = length(unique(x));
        un_exp = N * (1 - 1/e)^round;
        un_sim_list(round+1) = un_sim;
        un_exp_list(round+1) = un_exp;
    end
end
un_sim_list
un_exp_list

% Plots
figure();
scatter(rounds_list, un_sim_list);
hold on
scatter(rounds_list, un_exp_list);
legend("Simulation Results", "Theoretical Results");


% Auxiliary Functions
function y = BAGGING(x)
    LEN = length(x);
    for i = 1 : LEN
        y(i) = x(randi([1 LEN]));
    end
end



