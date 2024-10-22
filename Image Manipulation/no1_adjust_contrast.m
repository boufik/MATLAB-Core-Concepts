clc
clear all

% 1a. Image1 read
I = imread('littlegirl.png');
whos I
% 1b. Image1 show and histogram
figure(1);
imshow(I);
figure(2);
imhist(I);
display('Low contrast image: he intensity range of the image is rather narrow. The range does not cover the potential range of [0, 255], and is missing the high and low values that would result in good contrast.');

% 2a. Image2 creation - Contrast enhancement
% HISTOGRAM EQUALIZATION
I2 = histeq(I);
whos I2
% 2b. Image1 show and histogram
figure(3);
imshow(I2);
figure(4);
imhist(I2);
display('High contrast image');

