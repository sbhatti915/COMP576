% Sameer Bhatti
% sabhatti@live.unc.edu
% 2/10/20
% Assignment2.m
%
% Adds two images and displays the result

clc
clear
close all

%% Problem 1
a = imread('stripes.jpg');

a = im2double(a);

a = sqrt(a(:,:,1).^2 + a(:,:,2).^2 + a(:,:,3).^2);

IA = zeros(128, 128);
b = 0;

for x = 1:128
    for y = 1:128
        m = size(a);
        p = floor(m(1)/128)*x;
        z = floor(m(2)/128)*y;
        IA(x,y) = a(p,z);
        
    end
end


for j = 1:128
    for k = 1:128
        if k > j 
            IB(j,k) = 1000;
        end
        if k < j
            IB(j,k) = -1000;
        end
        if k == j
            IB(j,k) = 0;
        end
    end
end

IB = binomial_blurring(IB);

for x = 1:128
    for y = 1:128
        b = b + IB(x,y)^2;
    end
end

b = b/128^2;
std_dev = sqrt(b)*0.01;

% for j = 1:128
%     for k = 1:128
%         IB(j,k) = IB(j,k) + gnoise(std_dev);
%     end
% end

%% Problem 2
% [mag, phase]=AmpPhaseDFT(IB);
% 
% a = zeros(65,128);
% b = zeros(65,128);
% vx = zeros(65,128);
% vy = zeros(65,128);
% 
% for j = 1:128
%     for k = 1:65
%         if j <= 65
%             x = (2*pi*(j-1)/128)^2;
%         else
%             x = (2*pi*(j-129)/128)^2;
%         end
%         a(k,j) = x;
%     end
% end
% 
% for k = 1:128
%     for j= 1:65
%         b(j,k)=(2*pi*(j-1)/128)^2;
%     end
% end
% 
% a_inv = exp(-a);
% b_inv = exp(-b);
% 
% prod = mag .* b_inv .* a_inv;
% 
% for j = 1:128
%     for k = 1:65
%         if j<=65
%             vx(k,j) = (j-1)/128;
%         else
%             vx(k,j) = (j-129)/128;
%         end
%     end
% end
% 
% 
% for j = 1:128
%     for k = 1:65
%         vy(k,j) = (k-1)/128;
%     end
% end
% 
% vx(1,1) = 0;
% vx(65,1) = 0;
% vx(1,65) = 0;
% vx(65,65) = 0;
% 
% vy(1,1) = 0;
% vy(65,1) = 0;
% vy(1,65) = 0;
% vy(65,65) = 0;
% 
% Fx = prod .* vx;
% Fy = prod .* vy;
% 
% phase = phase + (pi/2);
% phase(1,1) = 0;
% phase(65,1) = 0;
% phase(1,65) = 0;
% phase(65,65) = 0;
% 
% IA = ReconfromAmpPhase(Fx,phase);
% IB = ReconfromAmpPhase(Fy,phase);

%% Problem 3
% [mag, phase]=AmpPhaseDFT(IB);
% 
% a = zeros(65,128);
% b = zeros(65,128);
% vx = zeros(65,128);
% vy = zeros(65,128);
% 
% for j = 1:65
%     for k = 1:128
%         if k <= 65
%             x = (2*pi*(k-1)/128)^2;
%         else
%             x = (2*pi*(k-129)/128)^2;
%         end
%         a(j,k) = x;
%     end
% end
% 
% for j = 1:65
%     for k = 1:128
%         b(j,k)=(2*pi*(j-1)/128)^2;
%     end
% end
% 
% a_inv = exp(-a);
% b_inv = exp(-b);
% 
% prod = mag .* b_inv .* a_inv;
% 
% for j = 1:128
%     for k = 1:65
%         if j<=65
%             vx(k,j) = (j-1)/128;
%         else
%             vx(k,j) = (j-129)/128;
%         end
%     end
% end
% 
% 
% for j = 1:128
%     for k = 1:65
%         vy(k,j) = (k-1)/128;
%     end
% end
% 
% vx(1,1) = 0;
% vx(65,1) = 0;
% vx(1,65) = 0;
% vx(65,65) = 0;
% 
% vy(1,1) = 0;
% vy(65,1) = 0;
% vy(1,65) = 0;
% vy(65,65) = 0;
% 
% Fx = prod .* vx;
% Fy = prod .* vy;
% 
% newPhase = phase + (pi/2);
% newPhase(1,1) = 0;
% newPhase(65,1) = 0;
% newPhase(1,65) = 0;
% newPhase(65,65) = 0;
% 
% Dx = ReconfromAmpPhase(Fx,newPhase);
% Dy = ReconfromAmpPhase(Fy,newPhase);
% 
% for j = 1:65
%     for k = 1:128
%         if k <=65
%             vxx(j,k) = -(((k-1)/128)^2);
%         else
%             vxx(j,k) = -(((k-129)/128)^2);
%         end
%     end
% end
% 
% for j = 1:65
%     for k = 1:128
%         vyy(j,k) = -(((j-1)/128)^2);
%     end
% end
% 
% vxy = -(vx .* vy);
% 
% vxx(1,1) = 0;
% vxx(1,65) = -0.25;
% vxx(65,65) = -0.25;
% vxx(65,1) = 0;
% 
% vyy(1,1) = 0;
% vyy(1,65) = 0;
% vyy(65,1) = -0.25;
% vyy(65,65) = -0.25;
% 
% vxy(1,1) = 0;
% vxy(1,65) = 0;
% vxy(65,1) = 0;
% vxy(65,65) = -0.25;
% 
% Fxx = prod .* vxx;
% Fxy = prod .* vxy;
% Fyy = prod .* vyy;
% 
% phase(1,1) = 0;
% phase(65,1) = 0;
% phase(1,65) = 0;
% phase(65,65) = 0;
% 
% xx = ReconfromAmpPhase(Fxx,phase);
% xy = ReconfromAmpPhase(Fxy,phase);
% yy = ReconfromAmpPhase(Fyy,phase);
% 
% for i = 1:128
%     for j = 1:128
%         pix = [xx(i,j), xy(i,j); xy(i,j), yy(i,j)];
%         [vec,val] = eig(pix);
%         eigenvalues(i,j,1) = val(1,1);
%         eigenvalues(i,j,2) = val(2,2);
%         gradient = [Dx(i,j),Dy(i,j)];
%         
%         if val(1,1) >= 0 && val(2,2) >= 0
%             vector = [];
%             image_e(i,j) = 0;
%         else
%             if val(1,1) < val(2,2)
%             vector = vec(:,1);
%             else
%             vector = vec(:,2);
%             end
%         image_e(i,j) = abs(dot(gradient,vector));
%         
%         end
%         vectors{i,j} = vector;
%     end
% end
% 
% image_e(:,1) = 0; 
% image_e(1,:) = 0; 
% image_e(128,:) = 0; 
% image_e(:,128) = 0;

%% Problem 4
[mag, phase]=AmpPhaseDFT(IB);

a = zeros(65,128);
b = zeros(65,128);
vx = zeros(65,128);
vy = zeros(65,128);
vector = zeros(128,128);

for j = 1:65
    for k = 1:128
        if k <= 65
            x = (2*pi*(k-1)/128)^2;
        else
            x = (2*pi*(k-129)/128)^2;
        end
        a(j,k) = x;
    end
end

for j = 1:65
    for k = 1:128
        b(j,k)=(2*pi*(j-1)/128)^2;
    end
end

a_inv = exp(-a);
b_inv = exp(-b);

prod = mag .* b_inv .* a_inv;

for j = 1:128
    for k = 1:65
        if j<=65
            vx(k,j) = (j-1)/128;
        else
            vx(k,j) = (j-129)/128;
        end
    end
end


for j = 1:128
    for k = 1:65
        vy(k,j) = (k-1)/128;
    end
end

vx(1,1) = 0;
vx(65,1) = 0;
vx(1,65) = 0;
vx(65,65) = 0;

vy(1,1) = 0;
vy(65,1) = 0;
vy(1,65) = 0;
vy(65,65) = 0;

Fx = prod .* vx;
Fy = prod .* vy;

newPhase = phase + (pi/2);
newPhase(1,1) = 0;
newPhase(65,1) = 0;
newPhase(1,65) = 0;
newPhase(65,65) = 0;

Dx = ReconfromAmpPhase(Fx,newPhase);
Dy = ReconfromAmpPhase(Fy,newPhase);

for j = 1:65
    for k = 1:128
        if k <=65
            vxx(j,k) = -(((k-1)/128)^2);
        else
            vxx(j,k) = -(((k-129)/128)^2);
        end
    end
end

for j = 1:65
    for k = 1:128
        vyy(j,k) = -(((j-1)/128)^2);
    end
end

vxy = -(vx .* vy);

vxx(1,1) = 0;
vxx(1,65) = -0.25;
vxx(65,65) = -0.25;
vxx(65,1) = 0;

vyy(1,1) = 0;
vyy(1,65) = 0;
vyy(65,1) = -0.25;
vyy(65,65) = -0.25;

vxy(1,1) = 0;
vxy(1,65) = 0;
vxy(65,1) = 0;
vxy(65,65) = -0.25;

Fxx = prod .* vxx;
Fxy = prod .* vxy;
Fyy = prod .* vyy;

phase(1,1) = 0;
phase(65,1) = 0;
phase(1,65) = 0;
phase(65,65) = 0;

xx = ReconfromAmpPhase(Fxx,phase);
xy = ReconfromAmpPhase(Fxy,phase);
yy = ReconfromAmpPhase(Fyy,phase);

for i = 1:128
    for j = 1:128
        pix = [xx(i,j), xy(i,j); xy(i,j), yy(i,j)];
        [vec,val] = eig(pix);
        eigenvalues(i,j,1) = val(1,1);
        eigenvalues(i,j,2) = val(2,2);
        gradient = [Dx(i,j),Dy(i,j)];
        
        if val(1,1) >= 0 && val(2,2) >= 0
            image_e(i,j) = 0;
        else
            if val(1,1) < val(2,2)
            vector = vec(:,1);
            else
            vector = vec(:,2);
            end
        image_e(i,j) = abs(dot(gradient,vector));
        
        end
        vectors(i,j,1) = vector(1);
        vectors(i,j,2) = vector(2);
    end
end

image_e(:,1) = 0; 
image_e(1,:) = 0; 
image_e(128,:) = 0; 
image_e(:,128) = 0;

for i = 2:127
    for j = 2:127
        y(i,j) = atan(vectors(i,j,2)/vectors(i,j,1));
        
        if y(i,j) >= pi/8 && y(i,j) <= 3*pi/8
            if image_e(i,j) > image_e(i-1,j+1) && image_e(i,j) > image_e(i+1,j-1)
                R_matrix(i,j) = 1;
            else
                R_matrix(i,j) = 0;
            end
            
        elseif y(i,j) < pi/8 && y(i,j) >= -pi/8
            ang = 0; % Horizontal
            if image_e(i,j) > image_e(i,j+1) && image_e(i,j) > image_e(i,j-1)
                R_matrix(i,j) = 1;
            else
                R_matrix(i,j) = 0;
            end
            
        elseif y(i,j) < -pi/8 && y(i,j) >= -3*pi/8
            ang = -pi/4;
            if image_e(i,j) > image_e(i-1,j-1) && image_e(i,j) > image_e(i+1,j+1)
                R_matrix(i,j) = 1;
            else
                R_matrix(i,j) = 0;
            end
        else
            if image_e(i,j) > image_e(i+1,j) && image_e(i,j) > image_e(i-1,j)
                R_matrix(i,j) = 1;
            else
                R_matrix(i,j) = 0;
            end
        end
    end
end
        
        
        