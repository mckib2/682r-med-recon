%% Homework 3
% Nicholas McKibben
% ECEn 682R
% 2018-03-08

clear;
close all;

%% (1) Field Map
load('resphantom2.mat'); clear ans;
n = 160;
osf = 2;
kosf = 5;

% Display the two reconstructions
[ im1os,~ ] = gridkb(d1,ks,wt,n,osf,kosf,'image');
[ im2os,~ ] = gridkb(d2,ks,wt,n,osf,kosf,'image');

% Trim down to size and flip
idy = round((.5*n*(osf-1)+1):(.5*n*(osf+1)));
idx = fliplr(idy);
im1 = im1os(idx,idy);
im2 = im2os(idx,idy);

figure(1);
subplot(1,3,1);
imshow(abs(im1),[]);
title('d1 recon');

subplot(1,3,2);
imshow(abs(im2),[]);
title('d2 recon');

% From the two data sets, compute a field map measured in Hz
fm = compute_fm(im1,te1,im2,te2);

% Limit our attention to pixels whose amplitude is greater than 10% of the
% maximum
msk = double(abs(im1) > .1*max(max(abs(im1))));

subplot(1,3,3);
imshow(abs(fm),[]);
title('FM (Hz)');
    
% What range of frequencies do we observe?
fprintf('Range of frequencies in field map: %f -> %f\n', ...
    min(min(fm.*msk)),max(max(fm.*msk)));

%% (2) Multifrequency Reconstruction
% Reconstruct the first data set for f equal to -128 Hz to 128 Hz in steps
% of 16 Hz.

fmin = -128;
fmax = 128;
fstep = 16;
tad = size(d1,1)*samp;
im_mf = mf_recon(d1,ks,wt,n,te1,tad,fmin,fmax,fstep);
fs = fmin:fstep:fmax;

% Plot at desired frequencies
des_fs = [ -128 -64 0 64 128 ];
figure(2);
for ii = 1:numel(des_fs)
    subplot(1,numel(des_fs),ii);
    imshow(abs(im_mf(:,:,fs == des_fs(ii))),[]);
end

%% (3) Field Map Based Reconstruction
fms = round((fm - fmin)/fstep) + 1;
fms = max(fms,1);
fms = min(fms,numel(fs));

im_mp = zeros(n,n);
for ii = 1:n
    for jj = 1:n
        im_mp(ii,jj) = im_mf(ii,jj,fms(ii,jj));
    end
end

figure(3);
imshow(abs(im_mp),[]);
title('Map Based Reconstruction');