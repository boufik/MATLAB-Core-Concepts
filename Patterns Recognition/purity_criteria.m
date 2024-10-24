clear all
clc

% Data
N = 1000;
c1_list = 0 : N;
c2_list = N - c1_list;
LEN = length(c1_list);
GINI_list = zeros(1, LEN);
entropy_list = zeros(1, LEN);
error_list = zeros(1, LEN);
p1_list = zeros(1, LEN);

for i = 1 : LEN
    c1 = c1_list(i);
    c2 = c2_list(i);
    p1 = c1 / (c1+c2);
    p2 = c2 / (c1+c2);
    y1 = GINI(p1, p2);
    y2 = entropy(p1, p2);
    y3 = error(p1, p2);
    GINI_list(i) = y1;
    entropy_list(i) = y2;
    error_list(i) = y3;
    p1_list(i) = p1;
end

figure();
plot(p1_list, GINI_list);
title("Comparison of metrics");
xlabel("p_1");
ylabel("Metrics' value");
hold on
plot(p1_list, entropy_list);
hold on
plot(p1_list, error_list);
hold on
legend("GINI", "Entropy", "Error");

% Auxiliary Functions - Metrics
function y = GINI(p1, p2)
    y = 1 - p1^2 - p2^2;
end

function y = entropy(p1, p2)
    y = - p1*log2(p1) - p2*log2(p2);
end

function y = error(p1, p2)
    y = 1 - max(p1, p2);
end