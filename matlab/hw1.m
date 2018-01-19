%% Homework 1
% Nicholas McKibben
% EE369C
% 2018-01-15

clear;
close all;

%% Introduction
d = [zeros(1,10),10:-1:1,0,1:10,zeros(1,10)];
x = -20:20;

figure(1);
subplot(2,1,1);
stem(x,d,'k');
xlabel('x');
ylabel('s(x)');

%% (1) Sinc Interpolation
% Take the test signal, upsample by factor of 10
xi = -20:0.1:20;
di = sinc_interp(d,x,xi);

subplot(2,1,2);
plot(xi,di,'k');
hold on;
stem(x,d,'k');
xlabel('x');
ylabel('s(x)');
title('Original Data and sinc interpolation');

%% (2) Bunched Samples
% Generate the bunched sample data by subsampling the sinc interpolated
% data
db = [di(3:20:end) di(8:20:end)];
xb = [xi(3:20:end) xi(8:20:end)];

figure(2);
subplot(2,1,1);
plot(xi,di,'k');
hold on;
stem(xb,db,'k')
xlabel('x');
ylabel('s(x)');
title('Bunched Samples');
legend('Original Sinc Interpolation','Bunched Samples');

% Now see if we can recover the original samples from the bunched samples
du = sinc_resample(db,xb,x);
subplot(2,1,2);
plot(x,du,'k');
hold on;
plot(x,d,'k--','LineWidth',2);
stem(xb,db,'k');
xlabel('x');
ylabel('s(x)');
title('Bunched Samples and sinc resampling');
legend('Recovered','Original','Bunched Samples');

%% (3) Random Samples
% Choose random samples from the sinc interpolated signal
ndx = randsample(1:length(di),50);
dr = di(ndx);
xr = xi(ndx);

figure(3);
subplot(2,1,1);
stem(xr,dr,'k');
title('Random Samples');

dur = sinc_resample(dr,xr,x);
subplot(2,1,2);
plot(x,dur,'k');
hold on;
plot(x,d,'k--','LineWidth',2);
stem(xr,dr,'k');
xlabel('x');
ylabel('s(x)');
title('Random Samples and sinc resampling');
legend('Recovered','Original','Random Samples');

%% (4) Interpolators

figure(4);
% Pick a couple of interesting random samples and plot the results
dp1 = zeros(1,length(dr));
dp1(floor(end/4)) = 1;
dp1u = sinc_resample(dp1,xr,x);

dp2 = zeros(1,length(dr));
dp2(floor(3*end/4)) = 1;
dp2u = sinc_resample(dp2,xr,x);

subplot(2,1,1);
stem(xr,ones(1,length(xr)));
hold on;
plot(xi,sinc_interp(dp1u.',x,xi));
plot(xi,sinc_interp(dp2u.',x,xi));
title('Random Sampling Interpolators');

% There are only two interpolators for the bunched samples case
dp1 = zeros(1,length(db));
dp1(end/4) = 1;
dp1u = sinc_resample(dp1,xb,x);

dp2 = zeros(1,length(db));
dp2(3*end/4) = 1;
dp2u = sinc_resample(dp2,xb,x);

subplot(2,1,2);
stem(xb,ones(1,length(xb)));
hold on;
plot(xi,sinc_interp(dp1u.',x,xi));
plot(xi,sinc_interp(dp2u.',x,xi));
title('Bunched Sampling Interpolators');

%% (5) Signal Noise and Sample Timing Jitter
% First consider noise (error in sampled values) and 50 random samples
drn = dr + 0.25*randn(size(dr)); % noise, sigma=.25
du = sinc_resample(drn,xr,x);

figure(5);
subplot(2,1,1);
plot(xi,di,'k'); % ground truth
hold on;
stem(x,du,'k'); % uniform samples
plot(xr,drn,'r.'); % random samples, with noise
title('50 Noisy Random Samples');
legend('Original','Uniform Resampled','Random Samples');
xlabel('This didn''t work so well...');

% Now try with 100 random samples
ndx = randsample(1:length(di),100);
dr = di(ndx);
xr = xi(ndx);

drn = dr + 0.25*randn(size(dr));
du = sinc_resample(drn,xr,x);

subplot(2,1,2);
plot(xi,di,'k');
hold on;
stem(x,du,'k');
plot(xr,drn,'r.');
title('100 Noisy Random Samples');
legend('Original','Uniform Resampled','Random Samples');
xlabel('Oversampling did a lot better!');

% Now let's consider the effects of uncertainty in the timing of the
% sampling (jitter)
xrn = xr + 0.05*randn(size(xr)); % sigma=0.05
du = sinc_resample(dr,xrn,x);

figure(6);
subplot(2,1,2);
plot(xi,di,'k'); % ground truth
hold on;
stem(x,du,'k'); % uniform samples
plot(xrn,dr,'r.'); % random samples, with jitter
title('Jitter, 100 Random Samples');
legend('Original','Uniform Resampled','Random Samples');
xlabel('Again, pretty terrible...');

% Back to 50 random samples
ndx = randsample(1:length(di),50);
dr = di(ndx);
xr = xi(ndx);

xrn = xr + 0.05*randn(size(xr)); % sigma=0.05
du = sinc_resample(dr,xrn,x);

subplot(2,1,1);
plot(xi,di,'k'); % ground truth
hold on;
stem(x,du,'k'); % uniform samples
plot(xrn,dr,'r.'); % random samples, with jitter
title('Jitter, 50 Random Samples');
legend('Original','Uniform Resampled','Random Samples');
xlabel('Oversampling did a lot better!');

function di = sinc_interp(d,x,xi)
    %
    % inputs
    % d -- uniformly sampled data points, spaced by 1
    % x -- uniform sample locations
    % xi -- locations to evaluation for the sinc interpolation
    % outputs
    % di -- sinc interpolated values at locations xi
    
    X = 1; % always spaced by 1
    s = sinc((xi - x.')/X);
    di = d*s;
end

function du = sinc_resample(dn,xn,xu)
    %
    % inputs
    % dn -- non-uniformly sampled data points
    % xn -- non-uniform sample locations
    % xu -- uniform sample points, spaced by 1
    % outputs
    % du-- uniformly sampled data
    
    X = 1; % always spaced by 1
    E = sinc((xn - xu.')/X);
    du = E.'\dn.';
end