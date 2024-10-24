clear all
clc


% *************************************************************************
% *************************************************************************
% Theory
% *************************************************************************
% *************************************************************************
% max{ P(C | A1, A2, A3, A4) }      <-- Bayes -->
% max{ P(A1, A2, A3, A4 | C) * P(C) }
% P(A1, A2, A3, A4 | C) = P(A1|C) * P(A2|C) * P(A3|C) * P(A4|C)
% Y = Yes, N = No, S = Sometimes, M = Mammals, NM = Non-Mammals



% *************************************************************************
% *************************************************************************
% Database
% *************************************************************************
% *************************************************************************
chars = ["Name", "Give Birth", "Can Fly", "Live in water", "Have legs", "Class"];
names = ["Human", "Python", "Salmon", "Whale", "Frog", "Komodo", "Bat", "Pigeon", "Cat", "Leopard Shark", "Turtle", "Penguin", "Porcupine", "Eel", "Salamander", "Gila Monster", "Platypus", "Owl", "Dolphin", "Eagle"]';
births = ["Y", "N", "N", "Y", "N", "N", "Y", "N", "Y", "Y", "N", "N", "Y", "N", "N", "N", "N", "N", "Y", "N"]';
flys = ["N", "N", "N", "N", "N", "N", "Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Y", "N", "Y"]';
waters = ["N", "N", "Y", "Y", "S", "N", "N", "N", "N", "Y", "S", "S", "N", "Y", "S", "N", "N", "N", "Y", "N"]';
legss = ["Y", "N", "N", "N", "Y", "Y", "Y", "Y", "Y", "N", "Y", "Y", "Y", "N", "Y", "Y", "Y", "Y", "N", "Y"]';
classes = ["M", "NM", "NM", "M", "NM", "NM", "M", "NM", "M", "NM", "NM", "NM", "M", "NM", "NM", "NM", "M", "NM", "M", "NM"]';
matrix = [names births flys waters legss classes]
ROWS = size(matrix, 1);
COLS = size(matrix, 2);

un_births = unique(births);
un_flys = unique(flys);
un_waters = unique(waters);
un_legss = unique(legss);

% Random input
name = "Random Animal";
birth = un_births(randi([1 length(un_births)]));
fly = un_flys(randi([1 length(un_flys)]));
water = un_waters(randi([1 length(un_waters)]));
legs = un_legss(randi([1 length(un_legss)]));
input = [name birth fly water legs];
% input = [name "Y" "N" "Y" "N"]
disp(mat2str(input));

% First, I have to count how many mammals and non-mammals are there
MAM_VEC = matrix(:, COLS);
indices_M = find(MAM_VEC == "M");
indices_NM = find(MAM_VEC == "NM");
num_M = length(indices_M);
num_NM = length(indices_NM);
prob_M = num_M / ROWS;
prob_NM = num_NM / ROWS;






% *************************************************************************
% *************************************************************************
% Algorithm
% *************************************************************************
% *************************************************************************


% Let class 1 be the correct class ("M" = "Mammals")
% I will check in the indices of "indices_M"
for i = 1 : num_M
    index = indices_M(i);
    MAM_BIRTHS_VEC(i) = births(index);
end
num1_M = length(find(MAM_BIRTHS_VEC == input(2)));

for i = 1 : num_M
    index = indices_M(i);
    MAM_FLYS_VEC(i) = flys(index);
end
num2_M = length(find(MAM_FLYS_VEC == input(3)));

for i = 1 : num_M
    index = indices_M(i);
    MAM_WATERS_VEC(i) = waters(index);
end
num3_M = length(find(MAM_WATERS_VEC == input(4)));

for i = 1 : num_M
    index = indices_M(i);
    MAM_LEGSS_VEC(i) =legss(index);
end
num4_M = length(find(MAM_LEGSS_VEC == input(5)));
P1 = PROBABILITY(ROWS, prob_M, num_M, num1_M, num2_M, num3_M, num4_M, "Mammal");


% Let class 2 be the correct class ("NM" = "Non-Mammals")
for i = 1 : num_NM
    index = indices_NM(i);
    NMAM_BIRTHS_VEC(i) = births(index);
end
num1_NM = length(find(NMAM_BIRTHS_VEC == input(2)));

for i = 1 : num_NM
    index = indices_NM(i);
    NMAM_FLYS_VEC(i) = flys(index);
end
num2_NM = length(find(NMAM_FLYS_VEC == input(3)));

for i = 1 : num_NM
    index = indices_NM(i);
    NMAM_WATERS_VEC(i) = waters(index);
end
num3_NM = length(find(NMAM_WATERS_VEC == input(4)));

for i = 1 : num_NM
    index = indices_NM(i);
    NMAM_LEGSS_VEC(i) =legss(index);
end
num4_NM = length(find(NMAM_LEGSS_VEC == input(5)));
P2 = PROBABILITY(ROWS, prob_NM, num_NM, num1_NM, num2_NM, num3_NM, num4_NM, "Non-Mammal");
compare(P1, P2, input);


% AUXILIARY FUNCTIONS
function prob = PROBABILITY(ROWS, prob_M, num_M, num1_M, num2_M, num3_M, num4_M, STRING)
    if(STRING == "Mammal")
        abbr = " M";
    else
        abbr = "NM";
    end
    prob1 = num1_M / num_M;
    prob2 = num2_M / num_M;
    prob3 = num3_M / num_M;
    prob4 = num4_M / num_M;
    disp("**************** " + STRING + " ****************");
    P_A_M = prob1 * prob2 * prob3 * prob4;
    disp("Pr(A|" + abbr + ") = " + num2str(num1_M) + "/" + num2str(num_M) + " * " + num2str(num2_M) + "/" + num2str(num_M) + " * " + num2str(num3_M) + "/" + num2str(num_M) + " * " + num2str(num4_M) + "/" + num2str(num_M) + " = ");
    disp("         = " + num2str(prob1) + " * " + num2str(prob2) + " * " + num2str(prob3) + " * " + num2str(prob4) + " = " + P_A_M);
    disp("Pr(" + abbr + ") = " + num2str(num_M) + "/" + num2str(ROWS) + " = " + num2str(prob_M));
    prob = P_A_M * prob_M;
    display(' ');
    disp("Pr(A|" + abbr + ")Pr(" + abbr + ") = " + num2str(prob));
    display(' ');
    display(' ');
end

function compare(P1, P2, input)
    if P1 >= P2
        disp(mat2str(input) + " = Mammal");
        disp("because: " + num2str(P1) + " >= " + num2str(P2));
    else
        disp(mat2str(input) + " = Non-Mammal");
        disp("because: " + num2str(P1) + " < " + num2str(P2));
    end
end


