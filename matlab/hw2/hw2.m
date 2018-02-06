%% Homework 2
% Nicholas McKibben
% ECEn 672R
% 2018-02-05

clear;
close all;
load('rt_spiral_03.mat');

%% (1) Simple Gridding Reconstruction
n = 128;
m = grid1(d,k,[],n,[],[]);
im = ifftshift(ifft2(m));
figure(1);
subplot(2,2,1);
imshow(abs(im),[]);
title('No Weights');

% Add weights
m2 = grid1(d,k,w,n,[],[]);
im2 = ifftshift(ifft2(m2));
subplot(2,2,2);
imshow(abs(im2),[]);
title('With Weights');

%% (2) Oversampled Gridding Reconstruction
o = 2; % oversampling factor
m3 = grid1(d,k,w,n,o,[]);
im3 = ifftshift(ifft2(m3));
subplot(2,2,3);
imshow(abs(im3(end/4:3*end/4-1,end/4:3*end/4-1)),[]);
title('2x k-space oversampling');

%% (3) Deapodization Correction
% x = linspace(1/min(real(k(:))),1/max(real(k(:))),n*o);
% y = linspace(1/min(imag(k(:))),1/max(imag(k(:))),n*o);
x = linspace(-2,2,n*o);
y = x;
dx = x(2) - x(1);
dy = y(2) - y(1);
c = sinc(x*dx).^2'*sinc(y*dy).^2;

im4 = im3./c;
subplot(2,2,4);
imshow(abs(im4(end/4:3*end/4-1,end/4:3*end/4-1)),[]);
title('Deapodization');

figure(2);
subplot(2,1,1);
im3_sl = im3(floor(3*end/5),end/4:3*end/4);
plot(abs(im3_sl),'-r');
hold on;
im4_sl = im4(floor(3*end/5),end/4:3*end/4);
plot(abs(im4_sl),'--b');
title('Slice through phantom');
legend('No deapodization','With deapodization');

subplot(2,1,2);
plot(abs(im3_sl - im4_sl));
title('Difference of the two');

%% (4) Improved Kernel
W = 5; % from Beatty paper
beta = pi*sqrt(((W/o)*(o - 1/2))^2 - 0.8);
m5 = grid1(d,k,w,n,o,[W,beta]);
im5 = ifftshift(ifft2(m5));

figure(3);
subplot(2,1,1);
kernel = kaiser(n*o,beta);
plot(x,kernel);

figure(4);
imshow(abs(im5),[]);