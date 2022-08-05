%% Adaptive Signal Processing
%% 2.1 The Least Mean Squares Algorithm
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 17/03/22

% Question 2.1b

clear
clc
close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

N = 1000; % no.samples
its = 100; % no. iterations
steps = [0.01, 0.05]; % step sizes
% we are 'shaping' the noise spectrum to simulate the signal x(n)
% hence, x(n) - a1x(n-1) - a2x(n-2) = noise
a1 = -0.1;
a2 = -0.8;
a0 = 1;
a = [a0, a1, a2];
sig = 0.25; % standard deviation of the noise
wgnoise = sqrt(sig).*randn(N,1);
x =  filter(1,a,wgnoise);
leak = 0;
figure
hold on
% first experiment just for one realisation
colors = {'b','r'};
for j = 1: length(steps)
    
    [xPred, ws, errs] =  LMS(x,steps(j),leak,length(a)-1); % excluding a0
    subplot(1,2,1)
    plot(10*log10(errs.^2),'color',colors{j},'LineWidth',1.5)
    hold on
end
title("\textbf{Squared Prediction Error}: One Realisation",'interpreter','latex','fontsize',16);
xlabel("Time Index" ,'interpreter','latex','fontsize',14);
ylabel("Squared Prediction Error (dB)",'interpreter','latex','fontsize',14);
legend('$\mu$ = 0.01','$\mu$=0.05','interpreter','latex')
grid on; grid minor;

empMisAd = [];

allerrs = zeros(N,its);
% repeat for 100 iterations
for j = 1: length(steps)
    
    for i = 1: its
        
    wgnoise = sqrt(sig).*randn(N,1);
    x =  filter(1,a,wgnoise);
        
    [xPred, ws, errs] =  LMS(x,steps(j),leak,length(a)-1); % excluding a0

    allerrs (:,i)= errs;
    end
    
    subplot(1,2,2)
    plot(10*log10(mean(allerrs.^2,2)),'color',colors{j},'LineWidth',1.5)
    hold on
    
    % Question 2.1 c
    % assuming that steady state is reached in both instances at approx t =
    % 600
    empMisAd = [empMisAd, (mean(mean(allerrs(600:end,:).^2,1),2)/0.25 -1)];
end

title("\textbf{Squared Prediction Error}: 100 Realisations",'interpreter','latex','fontsize',16);
xlabel("Time Index" ,'interpreter','latex','fontsize',14);
ylabel("Squared Prediction Error (dB)",'interpreter','latex','fontsize',14);
legend('$\mu$ = 0.01','$\mu$=0.05','interpreter','latex')
grid on; grid minor;
set(gcf,'color','w')


%% Question 2.1d

close all
N = 1000; % no.samples
its = 100; % no. iterations
steps = [0.01, 0.05]; % step sizes
% we are 'shaping' the noise spectrum to simulate the signal x(n)
% hence, x(n) - a1x(n-1) - a2x(n-2) = noise
a1 = -0.1;
a2 = -0.8;
a0 = 1;
a = [a0, a1, a2];
sig = 0.25; % standard deviation of the noise
wgnoise = sqrt(sig).*randn(N,1);
x =  filter(1,a,wgnoise);
leaks = [0];

empMisAd = [];
allCoeffs = zeros(length(a)-1,N,its);
% repeat for 100 iterations
figure
hold on
for j = 1: length(steps)
    for k = 1: length(leaks)
        
        for i = 1: its

        wgnoise = sqrt(sig).*randn(N,1);
        x =  filter(1,a,wgnoise);

        [xPred, ws, errs] =  LMS(x,steps(j),leaks(k),length(a)-1); % excluding a0

        allCoeffs(:,:,i)= ws;
        end
        
        
        subplot(1,2,j)
        h(1) = plot(mean(allCoeffs(1,:,:),3),'b','LineWidth',1.5);
        hold on
        if k == length(leaks)
            yline(0.1,'LineWidth',1.5)
        end

        h(2) = plot(mean(allCoeffs(2,:,:),3),'r','LineWidth',1.5);
        hold on
        if k == length(leaks)
            yline(0.8,'LineWidth',1.5)
        end
        title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{1}$ and $\\hat{a}_{2}$}: 100 Realisations, $\\mu$ = %0.2f',steps(j)),'interpreter','latex','fontsize',14);
        xlabel("Time Index" ,'interpreter','latex','fontsize',14);
        ylabel("Weight Values",'interpreter','latex','fontsize',14);
        legend([h(1),h(2)],'$\hat{a}_{1}$','$\hat{a}_{2}$','interpreter','latex','fontsize',15)
        grid on; grid minor;
        set(gcf,'color','w')


    end


end


%% Question 2.1f
close all
N = 1000; % no.samples
its = 100; % no. iterations
steps = [0.01, 0.05]; % step sizes
% we are 'shaping' the noise spectrum to simulate the signal x(n)
% hence, x(n) - a1x(n-1) - a2x(n-2) = noise
a1 = -0.1;
a2 = -0.8;
a0 = 1;
a = [a0, a1, a2];
sig = 0.25; % standard deviation of the noise
wgnoise = sqrt(sig).*randn(N,1);
x =  filter(1,a,wgnoise);
leaks = [0.1, 0.5, 0.8];
figure
hold on
% first experiment just for one realisation
colors = {[0 0 1],[1 0 0],[1 0 1]};

empMisAd = [];

allCoeffs = zeros(length(a)-1,N,its);
% repeat for 100 iterations
for j = 1: length(steps)
    figure
    for k = 1: length(leaks)
        
        for i = 1: its

        wgnoise = sqrt(sig).*randn(N,1);
        x =  filter(1,a,wgnoise);

        [xPred, ws, errs] =  LMS(x,steps(j),leaks(k),length(a)-1); % excluding a0

        allCoeffs(:,:,i)= ws;
        end
        
        
        subplot(1,2,1)
        plot(mean(allCoeffs(1,:,:),3),'color',colors{k},'LineWidth',1.5);
        hold on
        if k == length(leaks)
            yline(0.1,'LineWidth',1.5)
        end
        title("\textbf{Time Evolution of $a_{1}$}: 100 Realisations",'interpreter','latex','fontsize',14);
        xlabel("Time Index" ,'interpreter','latex','fontsize',14);
        ylabel("Leaky LMS Filter Weights Value",'interpreter','latex','fontsize',14);
        legend('$\mu$ = 0.05, $\gamma$ =0.2', '$\mu$ = 0.05, $\gamma$ =0.5','$\mu$ = 0.05, $\gamma$ =0.8','interpreter','latex')
        grid on; grid minor;
        
        
        subplot(1,2,2)
        plot(mean(allCoeffs(2,:,:),3),'color',colors{k},'LineWidth',1.5);
        hold on
        if k == length(leaks)
            yline(0.8,'LineWidth',1.5)
        end
        title("\textbf{Time Evolution of $a_{2}$}: 100 Realisations",'interpreter','latex','fontsize',14);
        xlabel("Time Index" ,'interpreter','latex','fontsize',14);
        ylabel("Leaky LMS Filter Weight Value",'interpreter','latex','fontsize',14);
        legend('$\mu$ = 0.05, $\gamma$ =0.2', '$\mu$ = 0.05, $\gamma$ =0.5','$\mu$ = 0.05, $\gamma$ =0.8','interpreter','latex')
        grid on; grid minor;
        set(gcf,'color','w')


    end


end


