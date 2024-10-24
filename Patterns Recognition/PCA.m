clear all
close all
clc

% Data - Step 0 (page 18 from Lecture 7)
display('************************************************');
display('************************************************');
display('        0. Initial-original Data');
x11 = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1]'
x22 = [2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9]'
n = length(x11);

% Step 1 - Minus mean
display('        1. Zero-mean data');
avg1 = mean(x11);
avg2 = mean(x22);
x1 = x11 - avg1
x2 = x22 - avg2
display(' ');

figure();
plot(x11, x22, 'rx');
title("10 spots - Initial Data");
% xlim([-2 4]);
% ylim([-2 4]);
hold on
plot(x1, x2, 'bx');
xlabel("x");
ylabel("y");
legend("Original Data", "Zero-mean Data");


% Steps 2, 3 - Covariance matrix, Eigenvalues, Eigenvectors
display('        2. Find the covariance matrix S');
S = cov(x1, x2)
display(' ');

eig_values = eig(S);
l1 = eig_values(1);
l2 = eig_values(2);
lambda = [l1 l2];
display('        3. Find the eigenvalues and eigenvectors');
display(' ');
disp("l_1 = " + num2str(l1));
disp("l_2 = " + num2str(l2));
[eig_vectors, no_need] = eig(S);
u1 = eig_vectors(:, 1);
u2 = -eig_vectors(:, 2);   % I put minus because I want the same results with the lecture (the u2 here and that in lecture have oppposite coordinates)
eig_vectors = [u1 u2];
disp("u_1 = " + mat2str(u1') + "^T");
disp("u_2 = " + mat2str(u2') + "^T");
display(' ');

figure();
plot(x1, x2, 'x');
title("Zero-mean Data with VERTICAL eigenvectors");
xlabel("x");
ylabel("y");
hold on
klisi1 = u1(1) / u1(2);
klisi2 = u2(1) / u2(2);
x = linspace(-2, 2.5, 101);
y1 = klisi1 * x;
y2 = klisi2 * x;
plot(x, y1, 'green');
hold on
plot(x, y2, 'blue');
legend("Zero-mean Data", "y = " + num2str(klisi1) + "*x", "y = " + num2str(klisi2) + "*x");



% Step 5 - Transpose of eigenvectors and sorting
RowFeatureVector = eig_vectors';
if l1 < l2
    % Need to change
    temp = RowFeatureVector(1, :);
    RowFeatureVector(1, :) = RowFeatureVector(2, :);
    RowFeatureVector(2, :) = temp;
end
RowFeatureVector;
RowZeroMeanData = [x1 x2]';
disp("        4. Calculate FinalData (they are separated by 'x' and 'y' axis)");
FinalData = (RowFeatureVector * RowZeroMeanData)'
display(' ');
display(' ');
% These data ("FinalData") are now separated by the 2 original axis "x" and
% "y" and not the eigenvectors, so this was the purpose of the
% transformation. WE CAN SAY THAT THE EIGENVECTORS ARE NOW THE AXIS!!!!

figure();
plot(FinalData(:, 1), FinalData(:, 2), 'x');
title("Data transformed with 2 eigenvectors (now eig = axis 'x', 'y')");
hold on
xx = linspace(-2, 2.5, 101);
yy = 0 * xx;
plot(xx, yy, 'blue');
hold on
yyy = linspace(-2, 2.5, 101);
xxx = 0 * yyy;
plot(xxx, yyy, 'green');






% Retrieve the information
% RowFeatureVector                                          2 x 2
% RowZeroMeanData = [x1 x2]';                               2 x 10
% FinalData = (RowFeatureVector * RowZeroMeanData)'         10 x 2

% Case 1 - Use 2 EIGENVECTORS
display('************************************************');
display('************************************************');
display('        Retrieve information');
display(' ');
InitialData = [x11 x22]
RowOriginalData = RowFeatureVector' * FinalData';       % (2x2) * (2x10) = 2 x 10
RowOriginalData = RowOriginalData';                     % 10 x 2
RowOriginalData(:, 1) = RowOriginalData(:, 1) + avg1;
RowOriginalData(:, 2) = RowOriginalData(:, 2) + avg2;
display('        When I use 2 eigenvectors (correct): ');
RowOriginalData
display(' ' );

% Case 2 - Use 1 EIGENVECTOR, THIS WITH LAMBDA = 1.28
% COMPRESS INFO - I WILL KEEP ONLY THE 1ST ROW = 1ST EIGENVECTOR
RowCompressedVector = RowFeatureVector(1, :);                   % 1 x 2
FinalCompressedData = FinalData(:, 1);                          % 10 x 1
RowCompressedData = RowCompressedVector' * FinalCompressedData';% (2x1) * (1x10) = 2 x 10
RowCompressedData = RowCompressedData';                         % 10 x 2
RowCompressedData(:, 1) = RowCompressedData(:, 1) + avg1;
RowCompressedData(:, 2) = RowCompressedData(:, 2) + avg2;
display('        When I use o 1 eigenvector (wrong): ');
RowCompressedData


% Calculate information loss
losses = zeros(n, 2);
for i = 1 : n
    for j = 1 : 2
        losses(i, j) = abs(RowOriginalData(i, j) - RowCompressedData(i, j)) / RowOriginalData(i, j);
    end
end

losses;
mean_losses = mean(losses);
loss1_perc = 100 * mean_losses(1);
loss2_perc = 100 * mean_losses(2);
display('        Information loss when using 1 eigenvector');
display(' ');
disp("Loss_1 = " + num2str(loss1_perc) + "%    Loss_2 = " + num2str(loss2_perc) + "%");
total_loss = mean(mean_losses);
disp("Total information loss = " + num2str(100 * total_loss) + "%");

figure();
plot(RowOriginalData(:, 1), RowOriginalData(:, 2), 'bx');
hold on
plot(RowCompressedData(:, 1), RowCompressedData(:, 2), 'rx');
title("Information Loss when using 1 eigenvectors");
xlabel("x");
ylabel("y");
legend("Using u_1, u_2", "Using only u_2");

