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
x = linspace(-1/2,1/2,n*o); y = x;
c = sinc(x).^2'*sinc(y).^2;

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
W = 5;
beta = pi*sqrt(((W/o)*(o - 1/2))^2 - 0.8);
im5 = grid1(d,k,w,n,o,[W,1e-6]);

figure(3);
subplot(2,1,1);
kernel = kaiser(n*o,beta);
plot(x,kernel);
title('Kasier-Bessel Kernel');

subplot(2,1,2);
plot(abs(im5(floor(3*end/5),end/4:3*end/4)));
title('Slice to show Apodization correction');

figure(4);
subplot(2,1,1);
imshow(abs(im5),[]);
title('Kasier-Bessel Kernel, Deappodization')

subplot(2,1,2);
imshow(abs(im5(2*n/4:3*2*n/4,n*2/4:3*n*2/4)),[]);
title('1x');