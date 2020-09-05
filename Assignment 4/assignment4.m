% Sameer Bhatti
% sabhatti@live.unc.edu
% 2/10/20
% Assignment2.m
%
% Adds two images and displays the result

clc
clear
close all

% I = grayscale image (n x n) - but we already provided the image in gray to reduce your work, so I is the input image
% L = incident light distribution
% R = object reflectance distribution

% n = log2(size(image,1));
% m = n-3;
image = imread('batman.png');

image = im2double(image);

image = sqrt(image(:,:,1).^2 + image(:,:,2).^2 + image(:,:,3).^2);

n = log2(size(image,1));
m = n-3;

boxLength = 2^m;
numBox = 8;

for x = 1:numBox
    for y = 1:numBox
        image_box(x,y,:,:) = image(((x - 1)*64 + 1): x*numBox^2 , ((y - 1)*64 + 1):y*numBox^2);
    end
end

for x = 1:numBox
    for y = 1:numBox
        control_val(x,y) = mean2(image_box(x,y,:,:));
    end
end


%% Part 2



k = [-1,3,-3,1;3,-6,3,0;-3,0,3,0;1,4,1,0]; % kernal

I = image(97:416,97:416);
s = [8 7 6 5 4 3 2 1];
xSpline = zeros(8,320);
ySpline = zeros(320,320);

% B spline values
for x = 1:8
    for y = 2:6
        for v = 1:64
            tVal = [((2*v - 1)/128)^3 , ((2*v - 1)/128)^2 , ((2*v - 1)/128) , 1];
            xSpline(x,(y-2)*64 + v) = 1/6*tVal*k*[control_val(x,y-1) ; control_val(x,y) ; control_val(x,y+1) ; control_val(x,y+2)];
        end
    end
end

% % v*(y-1)
for x = 1:320
    for y = 2:6
        for v = 1:64
            tVal = [((2*v - 1)/128)^3 , ((2*v - 1)/128)^2 , ((2*v - 1)/128) , 1];
            ySpline((y-2)*64 + v,x) = 1/6*tVal*k*[xSpline(y-1,x) ; xSpline(y,x) ; xSpline(y+1,x) ; xSpline(y+2,x)];
        end
    end
end

L = ySpline;

for x = 1:320
    for y = 1:320
        R(x,y) = I(x,y)/L(x,y);
    end
end