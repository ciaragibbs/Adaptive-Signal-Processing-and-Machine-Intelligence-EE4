%%  Real World Signals: Respiratory Sinus Arrhythmia from RR-Intervals
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 19/02/22
close all
clear 
clc
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
%% Question 1.5a

% Apply the standard periodogram as well as the averaged periodogram with
% different window lengths (e.g.50s 140s to obtain the PSD of RRI)

% Plott the PSDs of RRI data from three trials separately


% load("RAW-ECG.mat") % need to load also to get the sampling frequency
% 
% 
% Trial_1 = data(770:2.5e5);
% Trial_2 = data(2.564e5:4.99e5);
% Trial_3 = data(5.228e5:7.473e5);
% [xRRI1,fsRRI1]=ECG_to_RRI(Trial_1, fs);
% [xRRI2,fsRRI2]=ECG_to_RRI(Trial_2, fs);
% [xRRI3,fsRRI3]=ECG_to_RRI(Trial_3, fs);
% 
% % This tells you that the frequency of RRI1 RRI2 and RRI2 is 4 Hz!!!^^^

%%
load("RRI-DATA.mat")
xRRI1_new = detrend(xRRI1);
xRRI2_new = detrend(xRRI2);
xRRI3_new = detrend(xRRI3);

fsRRI = 4; % same for everything

% Bartlett's method will be used 

allRRI = {xRRI1_new; xRRI2_new; xRRI3_new};

figure
hold on
for i = 1:3
    
    RRI_oi = allRRI{i};
    RRI_len = length(RRI_oi);
    [stanPSD,stanf] = pwelch(RRI_oi, hamming(RRI_len), 0,1024*2, 4, 'onesided');
   
    [pSD150,f150] = pwelch(RRI_oi, hamming(150*fsRRI), 0,1024, 4, 'onesided');
    [pSD100,f100] = pwelch(RRI_oi, rectwin(100*fsRRI), 0,1024, 4, 'onesided');
    [pSD50,f50] = pwelch(RRI_oi, rectwin(50*fsRRI), 0,1024, 4, 'onesided');
    [pSD10,f10] = pwelch(RRI_oi, rectwin(10*fsRRI), 0,1024, 4, 'onesided');
    subplot(3,2,2*i-1)
    plot(stanf,10*log10(stanPSD),'b','LineWidth',1.5) 
    ax = gca;
    ax.FontSize = 12;
    xlabel('Frequency (Hz)')
    ylabel('PSD(dB)')
    title(sprintf('Standard Periodogram: RRI%d',i),'fontsize',15)
    
    grid on 
    grid minor
    subplot(3,2,2*i)
    plot(f150,10*log10(pSD150),'r','LineWidth',1.5)
    hold on
    plot(f100,10*log10(pSD100),'m','LineWidth',1.5)
    hold on
    plot(f50,10*log10(pSD50),'c','LineWidth',1.5)
    hold on
    plot(f10,10*log10(pSD10),'color',[0.5 0 0.5],'LineWidth',1.5)
    ax = gca;
    ax.FontSize = 12;
    
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)')
    title(sprintf('Averaged Periodogram: RRI%d',i),'fontsize',15)
    grid on
    grid minor
end 
legend('150 s window','100 s window','50 s window','10 s window')
set(gcf,'color','w')

%% Question 1.5c

% Plot the AR spectrum estimate for the RRI signals for the three trials
% To find the optimal AR model order, experiment with your model order
% until you observe a peak in the spectrum (approximately) correspondingto
% the theoretrical respiration rate.

figure
hold on
app = [8,10,25];
for i = 1:3
    
    RRI_oi = allRRI{i};
    RRI_len = length(RRI_oi);
    [stanPSD,stanf] = pwelch(RRI_oi, RRI_len, 0,1024, 4, 'onesided');
    subplot(1,3,i)
    plot(stanf,10*log10(stanPSD),'b','LineWidth',1.5)
    hold on
    % Now testing AR models
    for j = app(i)
        
        [h,w] = pyulear(RRI_oi, j, 1024, 4);
        plot(w,10*log10(h),'r','LineWidth',1.5)
        hold on
        
    end
    ax = gca;
    ax.FontSize = 14;
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)')
    title(sprintf('RRI%d',i),'fontsize',16)
    switch i
        case 1
            legend('Standard','AR(8)')
        case 2
            legend('Standard','AR(11)')
        case 3
            legend('Standard','AR(25)')
    end
    grid on
    grid minor
end 
sgtitle('Spectral Estimation with AR Modelling','fontsize',18)
set(gcf,'color','w')
