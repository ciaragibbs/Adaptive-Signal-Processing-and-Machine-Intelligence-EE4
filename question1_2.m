%% Classical and Modern Spectrum Estimation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 28/02/22
% Question 1.2 - Periodogram-based methods applied to real-world data
clear;
clc;
close all;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');


%% Sunspot data - 1.2a
load sunspot.dat

sunDat = sunspot(:,2);
% remove mean
sunDatMean = sunDat - mean(sunDat);
% detrend
sunDatDT = detrend(sunDat);
% log and then subtract mean
sunDatLog = log(sunDat + eps) - mean(log(sunDat + eps));

% calculate the respective periodograms

% the pre-processing will involve the use of different windows
% 1. A rectangular window
[PSDorig,worig] = periodogram(sunDat,ones(1,length(sunDatMean)),[],1);
[PSDMean,wMean] = periodogram(sunDatMean,ones(1,length(sunDatMean)),[],1);
[PSDDT,wDT] = periodogram(sunDatDT,ones(1,length(sunDatDT)),[],1);
[PSDLog,wLog] = periodogram(sunDatLog,ones(1,length(sunDatLog)),[],1);
% If the value of the PSD is very small at a given frequency, set to 1,
% such that the 10log10 (dB) is set to 0 for plotting
PSDMean(PSDMean < 0.001) = 1;
PSDDT(PSDDT < 0.001) = 1;
PSDLog(PSDLog < 0.001) = 1;


figure
subplot(1,2,1)
plot(worig,10*log10(PSDorig),'b', 'Linewidth', 1.5)
hold on
plot(wMean,10*log10(PSDMean),'r', 'Linewidth', 1.5)
hold on
plot(wDT,10*log10(PSDDT),'m',  'Linewidth', 1.5)
plot(wLog,10*log10(PSDLog),'c', 'Linewidth',1.5)
ax = gca;
ax.FontSize = 14;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
title('Periodogram Sunspot Data Series - Rectangular Window', 'Fontsize', 16)
legend('Raw Data','Mean Removal','Detrend','Log + Mean Removal','fontsize',13.5)
grid on
grid minor
%
% 2. A hamming window
[PSDorig,worig] = periodogram(sunDat,ones(1,length(sunDatMean)),[],1);
[PSDMean,wMean] = periodogram(sunDatMean,hamming(length(sunDatMean)),[],1);
[PSDDT,wDT] = periodogram(sunDatDT,hamming(length(sunDatDT)),[],1);
[PSDLog,wLog] = periodogram(sunDatLog,hamming(length(sunDatLog)),[],1);
PSDMean(PSDMean < 0.001) = 1;
PSDDT(PSDDT < 0.001) = 1;
PSDLog(PSDLog < 0.001) = 1;


subplot(1,2,2)
plot(worig,10*log10(PSDorig),'b', 'Linewidth', 1.5)
hold on
plot(wMean,10*log10(PSDMean),'r', 'Linewidth', 1.5)
hold on
plot(wDT,10*log10(PSDDT),'m', 'Linewidth', 1.5)
plot(wLog,10*log10(PSDLog),'c', 'Linewidth',1.5)
ax = gca;
ax.FontSize = 14;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
title('Periodogram Sunspot Data Series - Hamming Window', 'Fontsize', 16)
legend('Raw Data','Mean Removal','Detrend','Log + Mean Removal','fontsize',13.5)
grid on
grid minor
set(gcf,'color','w')


%% EEG Data - 1.2b

load('EEG_Data_Assignment1.mat');
sampNo = length(POz);
xRange = [11:20];
wins = [1,5,10];
winSamples = [1,5,10]*fs;
POz = POz - mean(POz);
% fs is provided in the file

% Compute the standard PSD
[PSDstandard,fstandard] = pwelch(POz, ones(1,sampNo), 0,fs*10, fs, 'onesided');

% Compute the average periodograms fro the different window lengths
% Now can use pwelch with a non sample length window size

% For a rectangular filter window

[PSD1, f1] = pwelch(POz, rectwin(winSamples(1)),0,fs*10, fs, 'onesided');
[PSD5, f5] = pwelch(POz, rectwin(winSamples(2)),0,fs*10, fs, 'onesided');
[PSD10, f10] = pwelch(POz, rectwin(winSamples(3)),0,fs*10, fs, 'onesided');  

%
figure
subplot(2,2,1)
plot(fstandard, 10*log10(PSDstandard),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('Orginal EEG Periodogram','fontsize',15)
grid on
grid minor

%
subplot(2,2,2)
plot(f1,10*log10(PSD1),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 1s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,3)
plot(f5,10*log10(PSD5), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 5s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,4)
plot(f10,10*log10(PSD10), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 10s Windowing','fontsize',15)
grid on
grid minor

sgtitle('Averaged Periodogram for EEG Data - Rectangular Window', 'Fontsize', 20)
set(gcf,'color','w')



%%

% For a Hamming filter window

[PSD1, f1] = pwelch(POz, hamming(winSamples(1)),0,fs*10, fs, 'onesided');
[PSD5, f5] = pwelch(POz, hamming(winSamples(2)),0,fs*10, fs, 'onesided');
[PSD10, f10] = pwelch(POz, hamming(winSamples(3)),0,fs*10, fs, 'onesided');  

%
figure
subplot(2,2,1)
plot(fstandard, 10*log10(PSDstandard),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('Orginal EEG Periodogram','fontsize',15)
grid on
grid minor

%
subplot(2,2,2)
plot(f1,10*log10(PSD1),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 1s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,3)
plot(f5,10*log10(PSD5), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 5s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,4)
plot(f10,10*log10(PSD10), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 10s Windowing','fontsize',15)
grid on
grid minor


sgtitle('Averaged Periodogram for EEG Data - Hamming Window', 'Fontsize', 20)
set(gcf,'color','w')

%% Direct comparison of the standard periodogram with the 10s window 

figure
plot(fstandard, 10*log10(PSDstandard),'Color','blue','LineWidth',2)
hold on
plot(f10,10*log10(PSD10), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 16;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
xlim([0,100])
ylim([-150, -90])
legend('Standard','Averaged - 10s Hamming Window')
title('Standard vs Averaged Periodogram - EEG Data', 'Fontsize', 20)
grid on
grid minor
set(gcf,'color','w')