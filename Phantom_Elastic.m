% Copyright 2017 Kyungtaek Jun
% Licensed under the MIT License
% https://www.mit.edu/~amini/LICENSE.md

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Virtual Alignment Method for Regular Elastic Sample with Phantom Image %%%
%%% PLOS ONE, 2018                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear all, clc;

%%% the position of fixed point in the polar coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_r = 100;         % Radius from the center of Reconstrucion
fp_phi = 45;        % Angle (degree)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initial set up
max_angle = 180;        % rotating angle
projections = 360;      % number of projection images
dtheta = max_angle/projections;
theta = 0:dtheta:max_angle-dtheta;

%%% make phantom image by 500*500 pixels
pnx = 300;      % size of phantom
P = phantom('Modified Shepp-Logan',pnx);

%%% padding 100 pixels around the image
Pad = 200;
inx = pnx + 2*Pad;
iny = inx;
IS1 = zeros(inx);     % Initial set up for Image Sample
for j = Pad + 1:iny - Pad       % put phantom on the center
    for i = Pad + 1:inx - Pad
        temp_val = P(i-Pad, j-Pad);
        if temp_val < 0    % remove negative values for making videos
            IS1(i,j) = 0;
        else
            IS1(i,j) = temp_val;
        end
    end
end
figure, imshow(IS1), title('Padded Phantom Image');
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Set Up second Phantom Image %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IS2 = zeros(inx);     % Initial set up for Image Sample

qPad = Pad/4;
for j = qPad+1:iny-qPad
    for i = qPad+1:inx-qPad
        IS2(i,j) = IS1(floor((i-qPad-1)/2)+Pad+1,floor((j-qPad-1)/2)+Pad+1)/4;
    end
end

figure, imshow(IS2), title('2nd Phantom Image');
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Make Sinograms %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st Sinogram
[isino1,xp] = radon(IS1,theta);
[sino_nx, sino_ny] = size(isino1);

figure, imagesc(theta,xp,isino1)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Sinogram of Image Sample 1 (Ideal Sinogram)');
set(gcf,'color','w');

% 2nd Sinogram
[isino2,xp] = radon(IS2,theta);
figure, imagesc(theta,xp,isino2)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Sinogram of Image Sample 2 (Ideal Sinogram)');
set(gcf,'color','w');

% 3rd Mixed sinogram
half_ny = round(sino_ny/2);
half_nx = floor(sino_nx/2);
for j = 1:sino_ny
    for i = 1:sino_nx
        if j <= half_ny
            msino(i,j) = isino1(i,j);
        else
            msino(i,j) = isino2(i,j);
        end
    end
end

figure, imagesc(theta,xp,msino)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Mixed Sinogram (Unideal Sinogram)');
set(gcf,'color','w');
        

% 4 & 5th ideal sinograms with different size
for j = 1:sino_ny
    for i = 1:sino_nx
        nsino_ori(i,j) = 0;
        nsino_db(i,j) = 0;
    end
end

for j = 1:sino_ny
    for i = 1:sino_nx
        if j <= half_ny
            nsino_ori(i,j) = msino(i,j);
        else
            if i <= half_nx
                nsino_ori(i,j) = (msino(2*i-1,j) + msino(2*i,j));
            end
        end
    end
end

figure, imagesc(theta,xp,nsino_ori)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Unideal Sinogram - Original Size');
set(gcf,'color','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% T function for original size %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:projections
    for i = 1:sino_nx
        isino_ori(i,j) = 0;        % Ideal Sinogram with Original Size
    end
end

deg_to_rad = pi/180;
Virtural_RA = round(sino_nx/2);     % Change Virtural_RA value -> change the size of Reconstruction ( required change of sino_nx value )
for j = 1:projections
    pjtd_CA_x = 0;
    total_MAC = 0;
    for i = 1:sino_nx
        pjtd_CA_x = pjtd_CA_x + i*nsino_ori(i,j);
        total_MAC = total_MAC + nsino_ori(i,j);
    end
    
    pjtd_CA_x = pjtd_CA_x/total_MAC;                      % Calculate CA (CM)
    T = fp_r*cos((j*dtheta-fp_phi)*deg_to_rad);       % T function with Radius and Angle
    pjtd_CA_x = pjtd_CA_x + T;
    
    ceil_CA_x = ceil(pjtd_CA_x);
    floor_CA_x = floor(pjtd_CA_x);
    
    right_ratio = pjtd_CA_x - floor_CA_x;
    left_ratio = ceil_CA_x - pjtd_CA_x;
    if left_ratio == 0
        left_ratio = 1;
    end

    move = Virtural_RA - floor_CA_x;
    
    for i = 2:sino_nx
        if i-1+move <= 1 | i-1+move >= sino_nx
            continue
        else
            isino_ori(i-1+move,j) = nsino_ori(i-1,j)*left_ratio + nsino_ori(i,j)*right_ratio;
        end
    end
end
        
figure, imagesc(theta,xp,isino_ori)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title(['Ideal Sinogram of Original Specimen - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');

Recon_ori = iradon(isino_ori,dtheta,inx);
figure, imshow(Recon_ori)
title(['Ideally Focused Reconstruction with Original Size - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');

q_nx = floor(sino_nx/6);

% 5th Ideal Sinogram
for j = 1:sino_ny
    for i = 1:sino_nx
        if j > half_ny
            nsino_db(i,j) = msino(i,j);
        else
            if i <= half_nx
                nsino_db(2*i - 1,j) = (msino(q_nx+i,j))/2;
                nsino_db(2*i,j) = nsino_db(2*i - 1,j);
            end
        end
    end
end
figure, imagesc(theta,xp,nsino_db)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title('Unideal Sinogram - Double Size');
set(gcf,'color','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% T function for doubled size %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:projections
    for i = 1:sino_nx
        isino_db(i,j) = 0;        % Ideal Sinogram with Expanded Size
    end
end
    
for j = 1:projections
    pjtd_CA_x = 0;
    total_MAC = 0;
    for i = 1:sino_nx
        pjtd_CA_x = pjtd_CA_x + i*nsino_db(i,j);
        total_MAC = total_MAC + nsino_db(i,j);
    end
    
    pjtd_CA_x = pjtd_CA_x/total_MAC;                      % Calculate CA (CM)
    T = fp_r*cos((j*dtheta-fp_phi)*deg_to_rad);       % T function with Radius and Angle
    pjtd_CA_x = pjtd_CA_x + T;
    
    ceil_CA_x = ceil(pjtd_CA_x);
    floor_CA_x = floor(pjtd_CA_x);
    
    right_ratio = pjtd_CA_x - floor_CA_x;
    left_ratio = ceil_CA_x - pjtd_CA_x;
    if left_ratio == 0
        left_ratio = 1;
    end

    move = Virtural_RA - floor_CA_x;
    
    for i = 2:sino_nx
        if i-1+move <= 1 | i-1+move >= sino_nx
            continue
        else
            isino_db(i-1+move,j) = nsino_db(i-1,j)*left_ratio + nsino_db(i,j)*right_ratio;
        end
    end
end

figure, imagesc(theta,xp,isino_db)
colormap(hot)
colorbar
xlabel('Parallel Rotation Angle - \theta (degrees)');
ylabel('Parallel Sensor Position - x\prime (pixels)');
title(['Ideal Sinogram with Exapnded Specimen - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');

Recon_db = iradon(isino_db,dtheta,inx + 2*fp_r);
figure, imshow(Recon_db)
title(['Ideally Focused Reconstruction with Double Size - Fixted Point is on T function: R = ', int2str(fp_r), ' pixels, Angle = ', int2str(fp_phi), '\circ']);
set(gcf,'color','w');
