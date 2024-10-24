clear all
close all
clc
 
 
% *************************************************************************
% *************************************************************************
 
% Data - Generate random 2-D points
Ns = [2, 3, 5] * 10;
means = [0.5, 1, 1.5];
stds = [0.3, 0.4, 0.5];
LEN = length(Ns);
N = sum(Ns);                                    % Total # of points
[x, y] = generatePoints(Ns, means, stds);
N = 500;
x = rand(1, N);
y = rand(1, N);
points = [x; y];
colors = ["green", "blue", "yellow", "magenta", "red"];
 
 
 
 
% *************************************************************************
% *************************************************************************
 
% Algorithm's Initialization
K = 3;
min_x = min(x);
max_x = max(x);
min_y = min(y);
max_y = max(y);
% x's = 1st row       y's = 2nd row
% Column 'i' contains the coordinates of i-th element/point
centroids_x = min_x + rand(1, K) * (max_x - min_x);
centroids_y = min_y + rand(1, K) * (max_y - min_y);
centroids = [centroids_x; centroids_y];
centroids_x = [points(1, 1) points(1, floor(N/2)) points(1, N)];
centroids_y = [points(2, 1) points(2, floor(N/2)) points(2, N)];
centroids = [centroids_x; centroids_y];
figure();
plot(x, y, 'bo');
title("Iteration 0 - Random centroids");
xlabel("x");
ylabel("y");
hold on
plot(centroids_x, centroids_y, 'r+', 'MarkerSize', 10);
 
 
 
% *************************************************************************
% *************************************************************************
 
% Generate 2 useful matrices: 'distances' = K x N, where K = # of centroids
% and N = # of points, so each column 'j' contains the distance of j-th
% point with each of the k centroids
distances = zeros(K, N);
% And 'groups' = 1 x N, where each element is 1, 2, ..., or 'K' and shows
% to which group-centroid the element belongs
groups = zeros(1, N);
 
% First, I have to calculate the distances matrix between the existing
% centroids and the given random points
for n = 1 : N
    point = points(:, n);
    for k = 1 : K
        centroid = centroids(:, k);
        distances(k, n) = distance(centroid, point);
    end
    % After the completion of n-th column, I can determine the index of
    % centroid this item belongs to
    column = distances(:, n);
    [MIN, index] = min(column);
    groups(n) = index;
end
display('***********************************************************');
round = 0;
disp("        Iteration " + num2str(round));
centroids
% disp("Centroids = " + mat2str(centroids));
display(' ');
display(' ');
groups;
 
 
 
 
 
 
% *************************************************************************
% *************************************************************************
 
% Iterations with for-loop
MAX_ROUND = 10;
for round = 1 : MAX_ROUND
    
    display('***********************************************************');
    disp("        Iteration " + num2str(round));
    round;
    points_around_centroid = zeros(1, K);
    coord_around_centroid = zeros(2, K);
    % This vector will contain the indeces of points that belong to centroid 'k'
    % I have to find out the new centroids based on the 'groups' vector
    for n = 1 : N
        points_around_centroid(groups(n)) = points_around_centroid(groups(n)) + 1;
        coord_around_centroid(1, groups(n)) = coord_around_centroid(1, groups(n)) + points(1, n);
        coord_around_centroid(2, groups(n)) = coord_around_centroid(2, groups(n)) + points(2, n);
    end
    points_around_centroid;
    coord_around_centroid;
    % Now, my 2 vectors contain for each 'k': the sum of coordinates of the
    % points around this k-th centroid and the number of these points
    for k = 1 : K
        centroids(1, k) = coord_around_centroid(1, k) / points_around_centroid(k);
        centroids(2, k) = coord_around_centroid(2, k) / points_around_centroid(k);
    end
    centroids;
    % Maybe in some iterations a centroid has no corresponding points
    % around it, so:
    for k = 1 : K
        if isnan(centroids(1, k)) == 1 || isnan(centroids(2, k)) == 1
            disp("NaN for centroid " + num2str(k));
            centroids(1, k) = points(1, floor(N/2));
            centroids(2, k) = points(2, floor(N/2));
        end
    end
    centroids
    
    
    
    
   
    
    
    
    
    % Now, new centroids have determined - Update 'distances', 'groups'
    for n = 1 : N
        point = points(:, n);
        for k = 1 : K
            centroid = centroids(:, k);
            distances(k, n) = distance(centroid, point);
        end
        column = distances(:, n);
        [MIN, index] = min(column);
        groups(n) = index;
    end
    
    % disp("Centroids = " + mat2str(centroids));
    display('***********************************************************');
    display(' ');
    display(' ');
    
    
    
    
    
    figure();
    for n = 1 : N
        cluster_k = groups(n);
        plot(points(1, n), points(2, n), 'o', 'Color', colors(cluster_k));
        hold on
    end
    plot(centroids(1, :), centroids(2, :),  'r+', 'MarkerSize', 10);
    
end
 
 
 
 
 
 
 
 
% *************************************************************************
% *************************************************************************
 
% Auxiliary Functions
function plot_spots(x, y, TITLE)
    figure();
    plot(x, y, 'bo');
    title(TITLE);
    xlabel("x");
    ylabel("y");
end
 
function [x, y] = generatePoints(Ns, means, stds)
    if length(Ns) ~= length(means) || length(Ns) ~= length(stds) || length(means) ~= length(stds)
        x = linspace(0, 3, 4);
        y = linspace(0, 3, 4);
        display('Error');
    end
    x = [];
    y = [];
    for i = 1 : length(Ns)
        temp_x = rand(1, Ns(i)) * stds(i) + means(i);
        temp_y = rand(1, Ns(i)) * stds(i) + means(i);
        x = [x temp_x];
        y = [y temp_y];
    end
end
 
 
function d = distance(p1, p2)
    % p1 = [2 3]'    p2 = [1 4]'
    d = sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2);
end