% -----------------------------------------------------------------------%
clear all;
close all;
clc

%% Make pixel array
len = 512; %number of pixels in the array
cen = len/2 + 1;
dx = 5.0e-6;    %pixel spacing, m  
df = 1/(len*dx);   %sample spacing in the spatial frequency domain, cycles/m

%% Make input object
%photo_source = imread('poles2.jpg');
%object_photo = photo_source;

object = imread('ImA.jpg');
bin_object = rgb2gray(object);
xaxis = ((-len/2):(len/2-1))*dx;
yaxis = -xaxis;

%% make filter
fxaxis = ((-len/2):(len/2-1))*df;
fyaxis = -fxaxis;
[FX,FY] = meshgrid(fxaxis,fyaxis);  %2-D arrays hold fx location and fy location of all points
freq_rad = sqrt(FX.^2 + FY.^2);
maxfreq = (len/2-1)*df;

cutoff_freq1 = 0.1*maxfreq;
filter1 = double(freq_rad <= cutoff_freq1); % pinhole

cutoff_freq2 = 0.12*maxfreq;
filter2 = double(freq_rad >= cutoff_freq2);

cutoff_freq3 = 0.2*maxfreq;
filter3 = double(freq_rad <= cutoff_freq3);

filter4 = filter2.*filter3; % 3rd Fresnel zone

%% Horizontal Single Slit Aperture
h_single_slit = zeros(len,len);
h_halfwidth   = 80;   % pixels
h_halfheight  = 12;    % pixels
h_single_slit((cen-h_halfheight):(cen+h_halfheight),(cen-h_halfwidth):(cen+h_halfwidth)) = ...
                                          ones(2*h_halfheight+1,2*h_halfwidth+1);                      

%% Vertical Single Slit Aperture
v_single_slit = zeros(len,len);
v_halfwidth   = 12;   % pixels
v_halfheight  = 80;    % pixels
v_single_slit((cen-v_halfheight):(cen+v_halfheight),(cen-v_halfwidth):(cen+v_halfwidth)) = ...
                                          ones(2*v_halfheight+1,2*v_halfwidth+1);
                                      
%% Vertical Double Slit Aperture
v_double_slit = zeros(len,len);
v_halfwidth   = 12;   % pixels
v_halfheight  = 80;    % pixels
v_spacing     = 60;   % pixels
v_double_slit((cen-v_halfheight):(cen+v_halfheight),...
    ((cen-v_spacing/2)-v_halfwidth):((cen-v_spacing/2)+v_halfwidth)) = ...
                                          ones(2*v_halfheight+1,2*v_halfwidth+1);
                                      
v_double_slit((cen-v_halfheight):(cen+v_halfheight),...
    ((cen+v_spacing/2)-v_halfwidth):((cen+v_spacing/2)+v_halfwidth)) = ...
                                          ones(2*v_halfheight+1,2*v_halfwidth+1);

%% Horizontal Double Slit Aperture
h_double_slit = zeros(len,len);
h_halfwidth   = 80;   % pixels
h_halfheight  = 12;    % pixels
h_spacing     = 60;   % pixels
h_double_slit(((cen-h_spacing/2)-h_halfheight):((cen-h_spacing/2)+h_halfheight),...
    (cen-h_halfwidth):(cen+h_halfwidth)) = ...
                                          ones(2*h_halfheight+1,2*h_halfwidth+1);
                                      
h_double_slit(((cen+h_spacing/2)-h_halfheight):((cen+h_spacing/2)+h_halfheight),...
    (cen-h_halfwidth):(cen+h_halfwidth)) = ...
                                          ones(2*h_halfheight+1,2*h_halfwidth+1);

%% Fourier transform object
ftobj           = fftshift(fft2(fftshift(object(:,:,3))));
ft_single_slit  = fftshift(fft2(fftshift(h_single_slit)));
ft_double_slit  = fftshift(fft2(fftshift(h_double_slit)));

%% Filter implementation
ftimg1 = ftobj.*filter1; %Pinhole
ftimg2 = ftobj.*filter2;
ftimg4 = ftobj.*filter4; %3rd fresnel zone

ftimg_h_single_slit = ftobj.*h_single_slit;
ftimg_h_double_slit = ftobj.*h_double_slit;

ftimg_v_single_slit = ftobj.*v_single_slit;
ftimg_v_double_slit = ftobj.*v_double_slit;
%==========================================================
img1 = abs(fftshift(ifft2(fftshift(ftimg1))));% Pinhole
img4 = abs(fftshift(ifft2(fftshift(ftimg4))));% 3rd Fresnel Zone

img1_h_single_slit = abs(fftshift(ifft2(fftshift(ftimg_h_single_slit))));
img2_h_double_slit = abs(fftshift(ifft2(fftshift(ftimg_h_double_slit))));

img3_v_single_slit = abs(fftshift(ifft2(fftshift(ftimg_v_single_slit))));
img4_v_double_slit = abs(fftshift(ifft2(fftshift(ftimg_v_double_slit))));


%% Pinhole and 3rd Fresnel zone
figure('NumberTitle', 'off', 'Name', 'Pinhole and 3rd Fresnel zone');
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('parula');

subplot(2,4,1);
imagesc(xaxis,yaxis,filter1);axis('image');
xlabel('x, m');ylabel('y, m');
title('Pinhole');
subplot(2,4,2);
imagesc(fxaxis,fyaxis,img1);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('Pinhole - Image plane');
subplot(2,4,3);
mesh(fxaxis,fyaxis,img1); % Pinhole
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,4)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img1(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

subplot(2,4,5);
imagesc(xaxis,yaxis,filter4);axis('image');
xlabel('x, m');ylabel('y, m');
title('3rd Fresnel zone');
subplot(2,4,6);
imagesc(fxaxis,fyaxis,img4);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('3F zone - Image plane');
subplot(2,4,7);
mesh(fxaxis,fyaxis,img4); % 3rd Fresnel Zone
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,8)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img4(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

%% Horizontal Single and Double slits
figure('NumberTitle', 'off', 'Name', 'Horizontal Single and Double slits');
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('parula');

subplot(2,4,1);
imagesc(xaxis,yaxis,h_single_slit);axis('image');
xlabel('x, m');ylabel('y, m');
title('h_single_slit');
subplot(2,4,2);
imagesc(fxaxis,fyaxis,img1_h_single_slit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('img1_single_lit');
subplot(2,4,3);
mesh(fxaxis,fyaxis,img1_h_single_slit);
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,4)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img1_h_single_slit(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

subplot(2,4,5);
imagesc(xaxis,yaxis,h_double_slit);axis('image');
xlabel('x, m');ylabel('y, m');
title('h_double_slit');
subplot(2,4,6);
imagesc(fxaxis,fyaxis,img2_h_double_slit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('3F zone - Image plane');
subplot(2,4,7);
mesh(fxaxis,fyaxis,img2_h_double_slit);
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,8)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img2_h_double_slit(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

%% Vertical Single and Double slits
figure('NumberTitle', 'off', 'Name', 'Vertical Single and Double slits');
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('parula');

subplot(2,4,1);
imagesc(xaxis,yaxis,v_single_slit);axis('image');
xlabel('x, m');ylabel('y, m');
title('v single slit');
subplot(2,4,2);
imagesc(fxaxis,fyaxis,img3_v_single_slit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('img3 v single slit');
subplot(2,4,3);
mesh(fxaxis,fyaxis,img3_v_single_slit);
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,4)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img3_v_single_slit(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

subplot(2,4,5);
imagesc(xaxis,yaxis,v_double_slit);axis('image');
xlabel('x, m');ylabel('y, m');
title('v double slit');
subplot(2,4,6);
imagesc(fxaxis,fyaxis,img4_v_double_slit);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('v double slit');
subplot(2,4,7);
mesh(fxaxis,fyaxis,img4_v_double_slit);
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,4,8)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img4_v_double_slit(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');

%% Mesh
mesh = imread('mesh.jpg');
bin_mesh = rgb2gray(mesh);

bw_thresh_mesh = double(imbinarize(bin_mesh));
bw_thresh_mesh2 = uint8(255*mat2gray(bw_thresh_mesh));

ftimg_mesh = ftobj.*bw_thresh_mesh;
img5_mesh = abs(fftshift(ifft2(fftshift(ftimg_mesh)))); % Mesh

%%
figure('NumberTitle', 'off', 'Name', 'Mesh');
set(gcf, 'Units','Normalized','OuterPosition',[0 0 1 1]);
colormap('parula');

subplot(2,2,1)
imagesc(xaxis,yaxis,bw_thresh_mesh);axis('image');
xlabel('x, m');ylabel('y, m');
title('bw thresh mesh');
subplot(2,2,2)
imagesc(fxaxis,fyaxis,img5_mesh);axis('image');
xlabel('fx, cycles/m');ylabel('fy, cycles/m');
colorbar('EastOutside');
title('img5 mesh');
subplot(2,2,3)
meshc(fxaxis,fyaxis,img5_mesh);
xlabel('fx, cycles/m');ylabel('fy, cycles/m');zlabel('Intensity');
title('Intensity values');
subplot(2,2,4)
plot(xaxis,bin_object(cen,:));hold on;grid on;
plot(xaxis,img5_mesh(cen,:),'r');
legend('object','image');xlabel('x, m');ylabel('Intensity');
title('Slice through centers of object and image');