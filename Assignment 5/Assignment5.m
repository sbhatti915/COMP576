% Sameer Bhatti
% sabhatti@live.unc.edu
% 4/16/20
% Assignment5.m
%
% Adds two images and displays the result

clc
clear
close all

%% Problem 1
% c_x = 40;
% c_y = 40; 
% c_z = 40;
% 
% r_x = 3; 
% r_y = 15; 
% r_z = 7;
% img = zeros(128,128,128);
% 
% for x = 1:128
%     for y = 1:128
%         for z = 1:128
%             if (x + 0.5 - c_x)/r_x^2 + (y + 0.5 - c_y)/r_y^2 + (z + 0.5 - c_z)/r_z^2 - 1 > 0
%                 a = 0;
%             else
%                 a = 1;
%             end
%             
%             if (x + 0.5 - c_x)/r_x^2 + (y + 0.5 - c_y)/r_y^2 + (z - 0.5 - c_z)/r_z^2 - 1 > 0
%                 b = 0;
%             else
%                 b = 1;
%             end
%             
%             if (x + 0.5 - c_x)/r_x^2 + (y - 0.5 - c_y)/r_y^2 + (z + 0.5 - c_z)/r_z^2 - 1 > 0
%                 c = 0;
%             else
%                 c = 1;
%             end
%             
%             if (x + 0.5 - c_x)/r_x^2 + (y - 0.5 - c_y)/r_y^2 + (z - 0.5 - c_z)/r_z^2 - 1 > 0
%                 d = 0;
%             else
%                 d = 1;
%             end
%             
%             if (x - 0.5 - c_x)/r_x^2 + (y + 0.5 - c_y)/r_y^2 + (z + 0.5 - c_z)/r_z^2 - 1 > 0
%                 e = 0;
%             else
%                 e = 1;
%             end
%             
%             if (x - 0.5 - c_x)/r_x^2 + (y + 0.5 - c_y)/r_y^2 + (z - 0.5 - c_z)/r_z^2 - 1 > 0
%                 f = 0;
%             else
%                 f = 1;
%             end
%             
%             if (x - 0.5 - c_x)/r_x^2 + (y - 0.5 - c_y)/r_y^2 + (z - 0.5 - c_z)/r_z^2 - 1 > 0
%                 g = 0;
%             else
%                 g = 1;
%             end
%             
%             if (x - 0.5 - c_x)/r_x^2 + (y - 0.5 - c_y)/r_y^2 + (z + 0.5 - c_z)/r_z^2 - 1 > 0
%                 h = 0;
%             else
%                 h = 1;
%             end
%             
%             img(x,y,z) = a + b + c + d + e + f + g + h;
%         end
%     end
% end
% 
% % z_slice = img(:,:,z_index);
% 
% %% Part B
% g = [40,40,30]; 
% v = [0, 2/3,5^(1/2)/3]; 
% angle = 3*pi/2;
% 
% quaternion = [cos(angle/2), sin(angle/2).*v];
% 
% qinv = [cos(angle/2), -sin(angle/2).*v];
% 
% rotatedPos = zeros(128^3,3);
% 
% c = 0;
% 
% for x = 1:128
%     for y = 1:128
%         for z = 1:128
%             u = [0,x-g(1),y-g(2),z-g(3)];
%             q = quatmultiply(quatmultiply(quaternion,u),qinv);
%             
%             c = c+1;
%             
%             rotatedPos(c,1) = q(2) + g(1);
%             
%             rotatedPos(c,2) = q(3) + g(2);
%             
%             rotatedPos(c,3) = q(4) + g(3);
%         end
%     end
% end
% 
% %% Part c
% c_x = 8; 
% c_y = 8; 
% c_z = 8;
% r_x = 3;
% r_y = 4;
% r_z = 5;
% size = 16;
% I = getI(c_x, c_y, c_z, size); 
% axis = [1, 0, 0]; 
% theta = pi/4;
% 
% z_index = 6; 
% quaternion = [cos(theta/2), sin(theta/2).*axis];
% 
% qinv = [cos(theta/2), -sin(theta/2).*axis];
% 
% rotatedPos = zeros(size^3,3);
% 
% c = 0;
% 
% for x = 1:size
%     for y = 1:size
%         for z = 1:size
%             u = [0,x-c_x,y-c_y,z-c_z];
%             q = quatmultiply(quatmultiply(quaternion,u),qinv);
%             
%             c = c+1;
%             
%             xPrime(x,y,z) = q(2) + c_x;
%             
%             yPrime(x,y,z) = q(3) + c_y;
%             
%             zPrime(x,y,z) = q(4) + c_z;
%         end
%     end
% end
% 
% Iprime = interp3(I,xPrime,yPrime,zPrime);
% slice = Iprime(:,:,z_index);
% 
% 


%% Part d
% dx = 0;
% dy = 0;
% dz = 0;
% sumx = 0;
% 
% for x = 1:16
%     for y = 1:16
%         for z = 1:16
%             
%             dx = dx + I(x,y,z)*x;
%             dy = dy + I(x,y,z)*y;
%             dz = dz + I(x,y,z)*z;
%             
%             sumx = sumx + I(x,y,z);
%         end
%     end
% end
% 
% dx = dx/sumx;
% dy = dy/sumx;
% dz = dz/sumx;
% 
% 
% dxx = 0;
% dxy = 0;
% dyy = 0;
% dyz = 0;
% dzz = 0;
% dxz = 0;
% 
% sumx = 0;
% 
% for x = 1:16
%     for y = 1:16
%         for z = 1:16
%             
%             dxx = dxx + I(x,y,z)*(x-d(1))^2;
%             dxy = dxy + I(x,y,z)*(x-d(1))*(y-d(2));
%             dyy = dyy + I(x,y,z)*(y-d(2))^2;
%             dyz = dyz + I(x,y,z)*(y-d(2))*(z-d(3));
%             dzz = dzz + I(x,y,z)*(z-d(3))^2;
%             dxz = dxz + I(x,y,z)*(x-d(1))*(z-d(3));
%                       
%             sumx = sumx + I(x,y,z);
%         end
%     end
% end
% 
% dxx = dxx/sumx;
% dxy = dxy/sumx;
% dyy = dyy/sumx;
% dyz = dyz/sumx;
% dzz = dzz/sumx;
% dxz = dxz/sumx;
% 
% A = [dxx,dxy,dxz; dxy, dyy, dyz; dxz, dyz, dzz];
% 
% [vec,val] = eig(A);
% 
% quat = [vec',vec]; % 4x3
% 
% for x = 1:3
%     n = norm(quat(:,x));
% end
% 
% quat = quat./n;

%%
A = [3.2551    0.0000    0.0000 ;0.0000    3.9924    1.6851 ;0.0000    1.6851    3.9924];

[vec,val] = eig(A);

quat = vec;

quat(2:4,:) = quat;
quat(1,:) = 0;

ratio = [val(1,1)/3, val(2,2)/4, val(3,3)/5];