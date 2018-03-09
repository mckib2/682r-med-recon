%% Homework 3
% Nicholas McKibben
% ECEn 682R
% 2018-03-08

clear;
close all;

%% (1) Field Map
load('resphantom2.mat'); clear ans;

% Display the two reconstructions
[ im1,~ ] = gridkb(d1,ks,wt,160,2,5,'image');

% Do some flipping to get our head on straight
im1 = flipud(im1);

figure(1);
imshow(abs(im1),[]);