% Copyright 2017 Kyungtaek Jun
% Licensed under the MIT License
% https://www.mit.edu/~amini/LICENSE.md

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Virtual Alignment Method for Rigid Sample with Phantom Image %%%
%%% Scientific Reports, 2017                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear all, clc;

%%% the position of fixed point in the polar coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_r = 100;         % Radius from the center of Reconstrucion
fp_phi = 45;        % Angle (degree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initial set up
max_angle = 180;        % rotating angle
projections = 360;      % number of projection images
dtheta = max_angle/projections;
theta = 0:dtheta:max_angle-dtheta;

%%% make phantom image by 500*500 pixels
pnx = 500;      % size of phantom
P = phantom('Modified Shepp-Logan',pnx);

%%% padding 100 pixels around the image
Pad = 100;
inx = pnx + 2*Pad;
iny = inx;
IS = zeros(inx);     % Initial set up for Image Sample
for j = Pad + 1:iny - Pad       % put phantom on the center
    for i = Pad + 1:inx - Pad
        temp_val = P(i-Pad, j-Pad);
        if temp_val < 0    % remove negative values for making videos
            IS(i,j) = 0;
        else
            IS(i,j) = temp_val;
        end
    end
end
figure, imshow(IS), title('Padded Phantom Image');
set(gcf,'color','w');

%%% random errors
ER = Pad;       % range of errors
vt_er = randi([-ER ER],1,projections);      % vertial translation errors
pt_er = randi([-ER ER],1,projections);      % parallel translation erros

%%% set up videos
writerObj = VideoWriter('Phantom_During_Scanning.avi');
open(writerObj);
for k = 1:projections
    new_IS = imtranslate(IS,[vt_er(k), pt_er(k)]);
    writeVideo(writerObj, new_IS);
end
close(writerObj);

%%% set up non ideal sinogram
sino_nonideal = [];
for k = 1:projections
    new_IS = imtranslate(IS,[vt_er(k), pt_er(k)]);
    [sino_clmn,xp] = radon(new_IS,theta(k));
    sino_nonideal = [sino_nonideal sino_clmn];
end

figure, imagesc(theta,xp,sino_nonideal)
colormap(winter)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Non-ideal sinogram with translation errors');
set(gcf,'color','w');
    
%%% ideal sinogram by applying virtual alignment method
[snx, sny] = size(sino_nonideal);
for j = 1:sny
    for i = 1:snx
        sino_ideal(i,j) = 0;
    end
end
    
deg_to_rad = pi/180;
Virtural_RA = round(snx/2);    
for j = 1:projections
    pjtd_CA_x = 0;
    total_MAC = 0;
    for i = 1:snx
        pjtd_CA_x = pjtd_CA_x + i*sino_nonideal(i,j);
        total_MAC = total_MAC + sino_nonideal(i,j);
    end
    
    pjtd_CA_x = pjtd_CA_x/total_MAC;                      % Calculate CA (CM)
    T = fp_r*cos((j*dtheta-fp_phi)*deg_to_rad);       % T function with Radius and Angle
    pjtd_CA_x = pjtd_CA_x + T;
    
    ceil_pjtd_CA_x = ceil(pjtd_CA_x);
    floor_pjtd_CA_x = floor(pjtd_CA_x);
    right_ratio = pjtd_CA_x - floor_pjtd_CA_x;
    left_ratio = ceil_pjtd_CA_x - pjtd_CA_x;
    if left_ratio == 0
        left_ratio = 1;
    end
    
    move = Virtural_RA - floor_pjtd_CA_x;
    
    for i = 2:snx-1
        if i+1+move <= 1 | i+move >= snx
            continue
        else
            sino_ideal(i+move,j) = sino_nonideal(i,j)*left_ratio + sino_nonideal(i+1,j)*right_ratio;   %
        end
    end
end
        
figure, imagesc(theta,xp,sino_ideal)
colormap(winter)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title(['Ideal Sinogram - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');

output_size = inx;
Recon = iradon(sino_ideal,theta,output_size);
figure, imshow(Recon)
title(['Ideally Focused Reconstruction - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');

