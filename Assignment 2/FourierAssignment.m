% Sameer Bhatti
% sabhatti@live.unc.edu
% 2/10/20
% FourierAssignment.m
%
% Adds two images and displays the result

clc
clear
close all

%% Problem 1
a = imread('car.jpeg');
b = imread('street.jpg');

a = im2double(a);
b = im2double(b);

a = sqrt(a(:,:,1).^2 + a(:,:,2).^2 + a(:,:,3).^2);
b = sqrt(b(:,:,1).^2 + b(:,:,2).^2 + b(:,:,3).^2);

IA = zeros(128, 128);
IB = zeros(128, 128);

for x = 1:128
    for y = 1:128
        m = size(a);
        p = floor(m(1)/128)*x;
        z = floor(m(2)/128)*y;
        IA(x,y) = a(p,z);
        
        n = size(b);
        p = floor(n(1)/128)*x;
        z = floor(n(2)/128)*y;
        IB(x,y) = b(p,z);
        
    end
end

[magA,phaseA] = AmpPhaseDFT(IA);
[magB,phaseB] = AmpPhaseDFT(IB);
OA = ReconfromAmpPhase(magA,phaseA);
OB = ReconfromAmpPhase(magB,phaseB);

%% Problem 2
% a = imread('car.jpeg');
% b = imread('street.jpg');
% 
% a = im2double(a);
% b = im2double(b);
% 
% a = sqrt(a(:,:,1).^2 + a(:,:,2).^2 + a(:,:,3).^2);
% b = sqrt(b(:,:,1).^2 + b(:,:,2).^2 + b(:,:,3).^2);
% 
% IA = zeros(128, 128);
% IB = zeros(128, 128);
% 
% for x = 1:128
%     for y = 1:128
%         m = size(a);
%         p = floor(m(1)/128)*x;
%         z = floor(m(2)/128)*y;
%         IA(x,y) = a(p,z);
%         
%         n = size(b);
%         p = floor(n(1)/128)*x;
%         z = floor(n(2)/128)*y;
%         IB(x,y) = b(p,z);
%         
%     end
% end
% 
% [magA,phaseA] = AmpPhaseDFT(IA);
% [magB,phaseB] = AmpPhaseDFT(IB);
% 
% AmpAPhaseB = ReconfromAmpPhase(magA,phaseB);
% AmpBPhaseA = ReconfromAmpPhase(magB,phaseA);
% 
% TextAnswer = "Phase";
% 
% %% Problem 3
% a = imread('car.jpeg');
% 
% del_y = 30;
% del_x = 10;
% 
% a = im2double(a);
% 
% a = sqrt(a(:,:,1).^2 + a(:,:,2).^2 + a(:,:,3).^2);
% 
% IA = zeros(128, 128);
% AS = zeros(128, 128);
% AS1 = zeros(128, 128);
% 
% for x = 1:128
%     for y = 1:128
%         
%         m = size(a);
%         p = floor(m(1)/128)*x;
%         z = floor(m(2)/128)*y;
%         IA(x,y) = a(p,z);
%      
%     end
% end
% 
% [magA,phaseA] = AmpPhaseDFT(IA);
% 
% for k = 1:128
%     for l = 1:128
%        
%         i = mod(k+del_y-1,128)+1;
%         j = mod(l+del_x-1,128)+1;
%         
%         AS1(k,j) = IA(k,l);
%         AS(i,j) = IA(k,l);
%     end
% end
% 
% % AS = circshift(IA,[0 delX]);
% [magASx,phaseASx] = AmpPhaseDFT(AS1);
% x_mag_diff = magASx - magA;
% x_phase_diff = phaseASx - phaseA;
% 
% % AS = circshift(IA,[delY delX]);
% [magASxy,phaseASxy] = AmpPhaseDFT(AS);
% xy_mag_diff = magASxy - magA;
% magnitude = "Same Magnitude";
% xy_phase_diff = phaseASxy - phaseA;
% phase = "Different Phase";