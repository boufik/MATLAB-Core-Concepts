clc
clear all

% ******** Data ********
rgb = imread('circularobjects.png');
gray = rgb2gray(rgb);
sens1 = 0.95;
sens2 = 0.97;
polarity1 = 'dark';                 % Red circles
polarity2 = 'bright';               % Blue circles
y1 = simulate(rgb, gray, polarity1, sens1);
y2 = simulate(rgb, gray, polarity1, sens2);
y3 = simulate(rgb, gray, polarity2, sens1);
y4 = simulate(rgb, gray, polarity2, sens2);



function sensitivity = simulate(rgb, gray, pol, sens)
    
    display('*************************************');
    if strcmp(pol, 'dark')
        color = 'r';
    elseif strcmp(pol, 'bright')
        color = 'b';
    end
    polarity = pol
    sensitivity = sens
    % ******** Find circles ********
    % 1. 'imfindcircles' finds circular objects that are brighter than the background.
    % So, set the parameter 'ObjectPolarity' to 'dark' in imfindcircles to search for dark circles.
    % 2. imfindcircles has a parameter 'Sensitivity' which can be used to control this internal threshold,
    % and consequently, the sensitivity of the algorithm. A higher 'Sensitivity' value sets the detection
    % threshold lower and leads to detecting more circles.
    [centers1,radii1] = imfindcircles(rgb,[20 25],'ObjectPolarity', polarity, ...
        'Sensitivity', sens)
    [centers2,radii2] = imfindcircles(gray,[20 25],'ObjectPolarity', polarity, ...
        'Sensitivity', sens)


    % ******** Visualize circles ********
    figure;
    imshow(rgb);
    h1 = viscircles(centers1, radii1, 'Color', color);
    figure;
    imshow(gray);
    h2 = viscircles(centers2, radii2, 'Color', color);
    disp("Found " + num2str(length(centers1)) + " circles.");
    
end