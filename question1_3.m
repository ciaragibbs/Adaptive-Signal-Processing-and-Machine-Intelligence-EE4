%% Classical and Modern Spectrum Estimation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 28/02/22
% Question 1.3 - Periodogram-based methods applied to real-world data
clear;
clc;
close all;

%% Part a: biased and unbiased ACF estimates
% number of signal samples
sampNo = 1024;
WGN = wgn(1,sampNo,1);
sampVec = 0:sampNo-1;
fs = 0.5;
sineWave = sin(sampVec*0.1) +  WGN;
filtWGN = filter([1/4 1/4 1/4 1/4], 1, WGN);

% calculate biased  ACF estimates of the signal
[WGN_b, WGN_lag] = xcorr(WGN, 'biased');
[sineWave_b, sineWave_lag] = xcorr(sineWave, 'biased');
[filtWGN_b, filtWGN_lag] = xcorr(filtWGN, 'biased');

% calculate unbiased  ACF estimates of the signal
[WGN_ub, WGN_ulag] = xcorr(WGN, 'unbiased');
[sineWave_ub, sineWave_ulag] = xcorr(sineWave, 'unbiased');
[filtWGN_ub, filtWGN_ulag] = xcorr(filtWGN, 'unbiased');
% since this is filtered with a rectangular window - the Fourier transform
% of a rectangular window is a sinc function. 
% Sinc function sill result in negative values - lose the positive semi
% definitness 

% Compute the correlogram
% Since the lag component makes the data 0 centered - it needs to be
% reversed shifted using ifftshift
PSD_wgn_b = fftshift(real(fft(ifftshift(WGN_b))));
PSD_sine_b = fftshift(real(fft(ifftshift(sineWave_b))));
PSD_fwgn_b = fftshift(real(fft(ifftshift(filtWGN_b))));

PSD_wgn_ub = fftshift(real(fft(ifftshift(WGN_ub))));
PSD_sine_ub = fftshift(real(fft(ifftshift(sineWave_ub))));
PSD_fwgn_ub = fftshift(real(fft(ifftshift(filtWGN_ub))));

% Plot WGN
figure
subplot(3,2,1)
plot(WGN_lag, WGN_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(WGN_lag, WGN_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('ACF')
title('WGN ACF')
legend('Unbiased','Biased')
grid on
grid minor
subplot(3,2,2)
plot(WGN_lag, PSD_wgn_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(WGN_lag, PSD_wgn_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('PSD')
title('WGN Correlogram')
legend('Unbiased','Biased')
grid on
grid minor

% Plot Noisy Sine
subplot(3,2,3)
plot(sineWave_lag, sineWave_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(sineWave_lag, sineWave_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('ACF')
title('Noisy Sine ACF')
legend('Unbiased','Biased')
grid on
grid minor
subplot(3,2,4)
plot(sineWave_lag, PSD_sine_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(sineWave_lag, PSD_sine_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('PSD')
title('Noisy Sine Correlogram')
legend('Unbiased','Biased')
grid on
grid minor
%
subplot(3,2,5)
plot(filtWGN_lag, filtWGN_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(filtWGN_lag, filtWGN_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('ACF')
title('Filtered WGN ACF')
legend('Unbiased','Biased')
grid on
grid minor
subplot(3,2,6)
plot(filtWGN_lag, PSD_fwgn_ub, 'Color', 'blue','LineWidth',2)
hold on
plot(filtWGN_lag, PSD_fwgn_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-sampNo, sampNo])
ylabel('PSD ')
title('Filtered WGN Correlogram')
legend('Unbiased','Biased')
grid on
grid minor
set(gcf,'color','w')


%% Part b: Plotting the PSD 

sampFreq = 10;
N = 1024;
n = 0:1/(2*sampFreq):sampFreq;
test_signal = [1.5*sin(2*pi*n*1) + 0.8*sin(2*pi*n*0.3) + 0.6*sin(2*pi*n*1.5) zeros(1, N-length(n))];
all_psds = zeros(100, N*2-1);

figure
hold on
for i = 1:100

    corrupted_signal = test_signal + wgn(length(test_signal), 1, 1)';
    % Correlogram
    [all_its_b, all_its_lag] = xcorr(corrupted_signal, 'biased');
    all_psds(i,:) = fftshift(real(fft(ifftshift(all_its_b))));
    subplot(1,2,1)
    plot((all_its_lag/max(all_its_lag))*sampFreq, real(all_psds(i,:)),'color','c','LineWidth',1.5);
    hold on
end

% Mean
plot((all_its_lag/max(all_its_lag))*sampFreq, mean(real(all_psds)),'color','b','LineWidth',1.5);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency (Hz)')
ylabel('PSD')
set(gca,'fontsize',15)
title('Realisations and Mean of the Power Spectral Density','fontsize',16)
grid on 
grid minor
subplot(1,2,2)
% Standard Deviation
plot((all_its_lag/max(all_its_lag))*sampFreq, std(real(all_psds)),'color','r','LineWidth',1.5);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency (Hz)')
ylabel('PSD')
set(gca,'fontsize',15)
title('Standard Deviation of the Power Spectral Density','fontsize',16)
grid on
grid minor
set(gcf,'color','w')

%% Part c: Realisations in dB

sampFreq = 10;
N = 1024;
n = 0:1/(2*sampFreq):sampFreq;
test_signal = [1.5*sin(2*pi*n*1) + 0.8*sin(2*pi*n*0.3) + 0.6*sin(2*pi*n*1.5) zeros(1, N-length(n))];
all_psds = zeros(100, N*2-1);

figure
hold on
for i = 1:100

    corrupted_signal = test_signal + wgn(length(test_signal), 1, 1)';
    % Correlogram
    [all_its_b, all_its_lag] = xcorr(corrupted_signal, 'biased');
    all_psds(i,:) = fftshift(real(fft(ifftshift(all_its_b))));
    subplot(1,2,1)
    plot((all_its_lag/max(all_its_lag))*sampFreq, 10*log10(real(all_psds(i,:))),'color','c','LineWidth',1.5);
    hold on
end

% Mean
plot((all_its_lag/max(all_its_lag))*sampFreq, 10*log10(mean(real(all_psds))),'color','b','LineWidth',1.5);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency (Hz)')
ylabel('PSD(dB)')
ylim([-30,20])
title('Realisations and Mean of the Power Spectral Density','fontsize',16)
grid on
grid minor
subplot(1,2,2)
% Standard Deviation
plot((all_its_lag/max(all_its_lag))*sampFreq, std(10*log10(real(all_psds))),'color','r','LineWidth',1.5);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
ylim([1,8])
title('Standard Deviation of the Power Spectral Density','fontsize',16)
grid on
grid minor
set(gcf,'color','w')

%% Part d: Complex Exponential Realisations

N= 1024;
n1 = [20 30 40];
n2 = [45 50 55 100];
colors = {'b','r','m','c'};
figure
hold on
subplot(1,2,1)
for i = 1: length(n1)
    
    n = 0:n1(i);
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    corrupted_signal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    zero_padding = zeros(1, N-length(corrupted_signal));
    psd_data = abs((fft([corrupted_signal zero_padding]))./length(n));
    f_axis = [0:N-1]/N;
    plot(f_axis, 10*log10(psd_data),'color',colors{i},'LineWidth',1.5)
    hold on
end
ax = gca;
ax.FontSize = 14;
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB)')
title('PSD: Noise-Corrupted Complex Exponential Signal','fontsize',14)
xlim([0.2 0.5])
legend('N=20','N=35','N=40')
grid on
grid minor
subplot(1,2,2)
for i = 1: length(n2)
    
    n = 0:n2(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    corrupted_signal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    zero_padding = zeros(1, N-length(corrupted_signal));
    psd_data = abs((fft([corrupted_signal zero_padding]))./length(n));
    f_axis = [0:N-1]/N;
    plot(f_axis,10*log10(psd_data),'color',colors{i},'LineWidth',1.5)
    hold on
end
ax = gca;
ax.FontSize = 14;
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (dB)')
xlim([0.2 0.5])
title('PSD: Noise-Corrupted Complex Exponential Signal','fontsize',14)
legend('N=45','N=50','N=55','N=100')
grid on
grid minor
set(gcf,'color','w')

%% Part e: MUSIC Algorithm

close all
N = 256;
n = 0:(30-1);
test_signal = [ exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) zeros(1, N - length(n))];
p = 2;
all_psds = zeros(100, 256);
figure
hold on
subplot(1,2,1)
for i = 1:100

    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    corrupted_signal = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    [~, R] = corrmtx(corrupted_signal, 14, 'mod');
    [all_psds(i, :), F] = pmusic(R, p, [ ], 1);
    plot(F, all_psds(i, :), 'color','c','LineWidth',1.5)
    hold on;

end

% mean
plot(F, mean(all_psds),'color','b','LineWidth',1.5)
ax = gca;
ax.FontSize = 14;
ylabel('Power Spectral Density')
xlabel('Frequency ($\pi$ radians)','Interpreter','latex')
xlim([0.2 0.4])
title('MUSIC Pseudospectrum - Realisations and Mean','Interpreter','latex','fontsize',15)
grid on
grid minor
% standard deviation
subplot(1,2,2)
plot(F, std(all_psds),'color','r','LineWidth',1.5)
ax = gca;
ax.FontSize = 14;
ylabel('Power Spectral Density')
xlabel('Frequency ($\pi$ radians)','Interpreter','latex')
xlim([0.2 0.4])
title('MUSIC Psuedospectrum - Standard Deviation','Interpreter','latex','fontsize',15)
grid on
grid minor
set(gcf,'color','w')

