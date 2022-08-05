%% Adaptive Signal Processing
%% 2.3 Adaptive Noise Cancellation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 10/04/22
clear
close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

%% Question 2.3a

sampNo = 1000;
deltas = [1,2,3,4];
mu = 0.01; 
its = 100;
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))
% generating the clean component of s(n)
x = sin(0.01*pi.*(0:sampNo-1)); 

figure
hold on
for i = 1: length(deltas)
    i
    xhat_all = zeros(its,sampNo);
    s_all =  zeros(its,sampNo);
    MSPE = zeros(its,1);
    
    for j = 1: its
        j
        v = randn(1,sampNo); % has unit variance
        % generating the coloured noise component of s(n)
        eta = filter(b,a,v);
        % constructing the input signal to ALE LMS
        s = x + eta;
        s_all(j,:) = s;
        [xhat,w,err] = ALE_LMS(s,mu,0,M,deltas(i));
        xhat_all(j,:) = xhat;
        MPSE(j) = mean((x(100:end)-xhat(100:end)').^2);
        
    end
    subplot(1,4,i)
    p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
    hold on
    p2 = plot(xhat_all','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
    hold on
    p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
    MSPE4title = mean(MPSE);
    title(sprintf('MSPE = %0.3f, $\\Delta$ = %0.0f',MSPE4title,deltas(i)),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Sample Index (n)','fontsize',15)
    ylabel('(AU)','fontsize',15)
    set(groot,'defaultLegendInterpreter','latex');
    allps = [p1(1),p2(1),p3];
    legend(allps)
    grid on
    grid minor
    
end

sgtitle('Empirical Justification of Minimum Delay','fontsize',18)
set(gcf,'color','w')
    
%% Question 2.3b

Ms = [5,10,15,20];
deltas = 1:25;
colors = {'b','r','m','c'};
MSPE_od =  zeros(numel(Ms),numel(deltas));
figure
hold on
ax = gca;
ax.FontSize = 15;
set(groot,'defaultLegendInterpreter','latex');
for i = 1: numel(Ms) % iterating across filter orders
    
    for j = 1: numel(deltas) % iterating across delays
        j
        xhat_all = zeros(its,sampNo);
        s_all =  zeros(its,sampNo);
        MSPE = zeros(its,1);
        
        for k = 1: its % across 100 iterations to take an average
            
            v = randn(1,sampNo); % has unit variance
            % generating the coloured noise component of s(n)
            eta = filter(b,a,v);
            % constructing the input signal to ALE LMS
            s = x + eta;
            s_all(j,:) = s;
            [xhat,w,err] = ALE_LMS(s,mu,0,Ms(i),deltas(j));
            xhat_all(j,:) = xhat;
            MSPE(k) = mean((x(100:end)-xhat(100:end)').^2); % offset to not consider any transient due to convergence time  
        end
            MSPE_od(i,j) = mean(MSPE);
    end
    plot(MSPE_od(i,:),'Color',colors{i},'LineWidth',2)
    hold on
end
xlabel('Delay ($\\Delta$)','fontsize',15)
ylabel('MSPE','fontsize',15)
title('Relationship between Delay and MSPE','fontsize',18)
legend('M=5','M=10','M=15','M=20')
grid on
grid minor
set(gcf,'color','w')

% curve for delay versus MSPE
sampNo = 1000;
deltas = [1:25];
mu = 0.01; 
its = 100;
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))
% generating the clean component of s(n)
x = sin(0.01*pi.*(0:sampNo-1)); 
MSPE_delay = zeros(1,numel(deltas));

for i = 1: length(deltas)
    i
    xhat_all = zeros(its,sampNo);
    s_all =  zeros(its,sampNo);
    MSPE = zeros(its,1);
    
    for j = 1: its
        
        v = randn(1,sampNo); % has unit variance
        % generating the coloured noise component of s(n)
        eta = filter(b,a,v);
        % constructing the input signal to ALE LMS
        s = x + eta;
        s_all(j,:) = s;
        [xhat,w,err] = ALE_LMS(s,mu,0,M,deltas(i));
        xhat_all(j,:) = xhat;
        MSPE(j) = mean((x(100:end)-xhat(100:end)').^2);
        
    end
    MSPE_delay(i) =  mean(MSPE);
    
end

figure
hold on
plot(MSPE_delay,'Color',colors{1},'LineWidth',2)
hold on
ax = gca;
ax.FontSize = 15;
set(groot,'defaultLegendInterpreter','latex');
xlabel('Delay ($\Delta$)','fontsize',15)
ylabel('MSPE','fontsize',15)
title('Relationship between Delay and MSPE, M=5','fontsize',16)
grid on
grid minor
set(gcf,'color','w')

% curve for filter order versus MSPE
sampNo = 1000;
deltas = 3;
mu = 0.01; 
its = 100;
M = [1:20]; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))
% generating the clean component of s(n)
x = sin(0.01*pi.*(0:sampNo-1)); 
MSPE_delay = zeros(1,numel(M));

for i = 1: length(M)
    i
    xhat_all = zeros(its,sampNo);
    s_all =  zeros(its,sampNo);
    MSPE = zeros(its,1);
    
    for j = 1: its
        
        v = randn(1,sampNo); % has unit variance
        % generating the coloured noise component of s(n)
        eta = filter(b,a,v);
        % constructing the input signal to ALE LMS
        s = x + eta;
        s_all(j,:) = s;
        [xhat,w,err] = ALE_LMS(s,mu,0,M(i),deltas);
        xhat_all(j,:) = xhat;
        MSPE(j) = mean((x(100:end)-xhat(100:end)').^2);
        
    end
    MSPE_delay(i) =  mean(MSPE);
    
end

figure
hold on
plot(MSPE_delay,'Color',colors{2},'LineWidth',2)
hold on
ax = gca;
ax.FontSize = 15;
set(groot,'defaultLegendInterpreter','latex');
xlabel('Model Order (M)','fontsize',15)
ylabel('MSPE','fontsize',15)
title('Relationship between Model Order and MSPE, $\Delta$=3','fontsize',18)
grid on
grid minor
set(gcf,'color','w')

%% Question 2.3c

sampNo = 1000;
deltas = [3];% only one delay value to be tested
mu = 0.01; 
its = 100;
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))
% generating the clean component of s(n)
x = sin(0.01*pi.*(0:sampNo-1)); 

figure
hold on

s_all =  zeros(its,sampNo);
xhat_allALE = zeros(its,sampNo);
MSPE_ALE = zeros(its,1);

xhat_allANC = zeros(its,sampNo);
MSPE_ANC = zeros(its,1);

for j = 1: its
    j
    v = randn(1,sampNo); % has unit variance
    % generating the coloured noise component of s(n)
    eta = filter(b,a,v);
    % constructing the input signal to ALE LMS
    s = x + eta;
    u = 1.2*eta+0.1;
    s_all(j,:) = s;

    % ALE
    [xhatALE,wALE,errALE] = ALE_LMS(s,mu,0,M,deltas);
    xhat_allALE(j,:) = xhatALE;
    MSPE_ALE(j) = mean((x(100:end)-xhatALE(100:end)').^2);
    % ANC
    [noiseEst,wANC,xhatANC] = ANC_LMS(s,u,mu,0,M);
    xhat_allANC(j,:) = xhatANC;
    MSPE_ANC(j) = mean((x(100:end)-xhatANC(100:end)').^2);

end
% plot ALE
subplot(1,3,1)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allALE','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ALE);
title(sprintf('MSPE = %0.3f, $\\Delta$ = 3, M =5',MSPE4title),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on
grid minor 
% sgtitle('Empirical Justification of Minimum Delay','fontsize',18)
% plot ALE
subplot(1,3,2)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allANC','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ANC);
title(sprintf('MSPE = %0.3f, $\\Delta$ = 3, M =5',MSPE4title),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on
grid minor 
set(gcf,'color','w')

% ensemble of realisations
ANC_ensemble = mean(xhat_allANC);
ALE_ensemble = mean(xhat_allALE);
subplot(1,3,3)
plot(ALE_ensemble,'b','LineWidth',2);
hold on
plot(ANC_ensemble,'r','LineWidth',2);
hold on
plot(x,'y','LineWidth',2);
title('Ensemble Means: ALE vs ANC','Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample Index (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
legend('ALE','ALC','$x(n)$')
grid on
grid minor 
set(gcf,'color','w')

%% Question 2.3d
clear
load("EEG_Data_Assignment2.mat");
% create a time axis
t =  (0:1/fs:(1/fs)*(length(Cz)-1));
sineWave = sin(2*pi*50*t);
v = 0.005;
noise = sqrt(v)*randn(1,length(Cz));
ref = sineWave + noise;
% define learning rates and model orders to be tested
mus = [0.001,0.01,0.1];
Ms = [5,10,20];

windowLength = 4096;
overlap = 2/3;
nfft = 5*windowLength;

figure
spectrogram(Cz, hanning(windowLength), round(overlap * windowLength), nfft, fs, 'yaxis');
title('EEG: Cz Noise-Corrupted Spectrogram','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Time (mins)','fontsize',15)
ylabel('Frequency (Hz)','fontsize',15)
yticks(0:10:60);
ylim([0 60]);
c = colorbar('TickLabelInterpreter', 'latex');
c.Label.String = "Power (dB)";
c.Label.Interpreter = 'latex';
set(gcf,'color','w');

% iterating over learning rates and model orders, using ANC-LMS
count = 1;
figure
hold on
for i=1:numel(Ms)
    i
    for j=1:numel(mus)
        % ANC-LMS
        [noiseEst,wANC,xhatANC] = ANC_LMS(Cz,ref,mus(j),0,Ms(i));
%         xhat_allANC(j,:) = xhatANC;
%         MSPE_ANC(j) = mean((x(100:end)-xhatANC(100:end)').^2);    
         
        subplot(3,3,count)
        spectrogram(xhatANC, hanning(windowLength), round(overlap * windowLength), nfft, fs, 'yaxis');
        title(sprintf('Cz Spectrogram: $\\mathbf{M}$ = %0.0f and $\\mathbf{\\mu}$ = %0.3f',Ms(i),mus(j)),'Interpreter','latex','fontsize',12)
        ax = gca;
        ax.FontSize = 13;
        xlabel('Time (mins)','fontsize',13)
        ylabel('Frequency (Hz)','fontsize',13)
        yticks(0:10:60);
        ylim([0 60]);
        c = colorbar('TickLabelInterpreter', 'latex');
        c.Label.String = "Power (dB)";
        c.Label.Interpreter = 'latex';
        set(gcf,'color','w');
        count = count + 1;
    end
end

%% Periodogram for the optimum model order and learning rate
M = 10;
mu = 0.001;
%Computing ANC algorithm
[noiseEst,wANC,xhatANC] = ANC_LMS(Cz,ref,mu,0,M);
% define parameters
windowLength = 1024;
nfft = 5*windowLength;
% periodogram of Cz
[p_corrupt,w_corrupt] = pwelch(Cz,rectwin(windowLength),0,nfft,fs,'onesided');
p_corrupt = 10*log10(p_corrupt);
% periodogram of denoised Cz
[p_denoise,w_denoise] = pwelch(xhatANC,rectwin(windowLength),0,nfft,fs,'onesided');
p_denoise = 10*log10(p_denoise);

figure
subplot(1,2,1)
plot(w_corrupt,p_corrupt,'b','Linewidth',1.5)
hold on
plot(w_denoise,p_denoise,'r','Linewidth',1.5)
title('Cz: Periodogram Before and After Denoising','Fontsize',14)
ax = gca;
ax.FontSize = 13;
xlabel('Frequency (Hz)','fontsize',14)
xlim([0 70])
ylabel('Power Density(dB)','fontsize',14)
legend('Noise Corrupted','Denoised')
grid on
grid minor
subplot(1,2,2)
[p_corrupt,w_corrupt] = pwelch(Cz,rectwin(windowLength),0,nfft,fs,'onesided');
[p_denoise,w_denoise] = pwelch(xhatANC,rectwin(windowLength),0,nfft,fs,'onesided');
plot(w_denoise, 10*log10(abs(p_corrupt-p_denoise)),'m','LineWidth',1.5)
title('Cz: Error Periodogram by Denoising','Fontsize',14)
ax = gca;
ax.FontSize = 13;
xlabel('Frequency (Hz)','fontsize',14)
xlim([0 70])
ylabel('Power Density (dB)','fontsize',14)
grid on
grid minor
set(gcf,'color','w')
