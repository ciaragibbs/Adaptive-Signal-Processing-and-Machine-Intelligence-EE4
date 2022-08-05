%% Classical and Modern Spectrum Estimation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 28/02/22
% Question 1.1
clear all;
close all;
clc;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
%%

sampNo = 1024*3;
sampFreq = 1000;
sampPeriod = 1/sampFreq;
sampleVector = 0:sampNo-1;
%
pulseGen = zeros(1,sampNo);
pulseGen(sampNo/2) = 50;


sineGen = sin(2*pi*0.5*sampleVector*sampPeriod);
% 
acfPulse = xcorr(pulseGen,'biased');
acfSine = xcorr(sineGen,'biased');
% 

% Use definition 1: DTFT of ACF
psdSine1 = abs(fftshift(fft(acfSine)));
psdPulse1 = abs(fftshift(fft(acfPulse)));

% Use definition 2: 
psdSine2 = abs(fftshift(fft(sineGen))).^2 / sampNo;
psdPulse2 = abs(fftshift(fft(pulseGen))).^2 / sampNo;

% Build the frequency axis in radians
freqRad1 = [-pi+pi/length(acfPulse):2*pi/length(acfPulse):pi-pi/length(acfPulse)];
% Convert the frequency axis in Hz
freqHz1 = freqRad1./(2*pi).*sampFreq;


% Build the frequency axis in radians
freqRad2 = [-pi+pi/length(psdPulse2):2*pi/length(psdPulse2):pi-pi/length(psdPulse2)];
% Convert the frequency axis in Hz
freqHz2 = freqRad2./(2*pi).*sampFreq;

%%

figure;
subplot(3, 1, 1);
plot([-sampNo+1:sampNo-1],acfPulse,'b','LineWidth', 2);
hold on;
plot([-sampNo+1:sampNo-1],acfSine,'r', 'LineWidth', 2);
legend('Pulse','Sine');
ax = gca;
ax.FontSize = 15;
title('Signal ACFs','fontsize',15)
xlabel('Time (s)')
grid on
grid minor
% title('Trend of ACF');
xlabel('Lags (sample)');
ylabel('ACF');
subplot(3, 1, 2);
plot(freqHz1, psdSine1,'b', 'LineWidth', 2);
hold on;
plot(freqHz2, psdSine2,'r--', 'LineWidth', 2);
legend('Definition 1', 'Definition 2');
title('Periodogram of sinusoids','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Normalised frequency ($\pi$ rads/sample)');
ylabel('PSD');
grid on 
grid minor
subplot(3, 1, 3);
plot(freqHz1, psdPulse1,'b', 'LineWidth', 2);
hold on;
plot(freqHz2, psdPulse2, 'r--', 'LineWidth', 2);
legend('Definition 1', 'Definition 2');
title('Periodogram of impulse','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Normalised frequency ($\pi$ rads/sample)');
ylabel('PSD');
grid on
grid minor
set(gcf,'color','w')
