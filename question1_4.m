%% Estimation of Autoregressive Processes
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 12/04/22
close all
clear
clc
rng default
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
%%
% Question 1.b
n = 1000;
rem = 500;
a = [1, -2.76, 3.81, -2.65, 0.92];
x =  filter(1,a, randn(n,1));
% discard first 500 samples to remove transient filter effects
x = x(rem+1:end);
% obtain the TRUE power spectrum 

% freqz returns th n-point freq response, as well as the angular freq for a
% digitial filter defined by the function coefficients stored in b and a
% [h,w] = freqz(b,a,n) 
figure
hold on
[h,w] = freqz(1,a,length(x));
thePSD = abs(h).^2;
ordersOI = [2,4,11];
colors = {'r','m','c',[0.5,0,0.5],[0.5,0.5,0.5]};
MSEs0 = [];
for i = 1: length(ordersOI)
    
    % The Yule-Walker method provides an autoregressive power spectral
    % density estimate
    [aest,e] = aryule(x,ordersOI(i));
    [h,w] = freqz(e^(1/2),aest,length(x));
    estPSD = abs(h).^2;
    subplot(2,3,i)
    plot(w/pi,10*log10(thePSD),'b','LineWidth',2);
    hold on
    plot(w/pi, 10*log10(estPSD),'r','LineWidth',2)
    MSEs0 = [MSEs0,mean((10*log10(thePSD) - 10*log10(estPSD)).^2)];
    ax = gca;
    ax.FontSize = 14;
    axis tight
    xlabel('Normalised Frequency ($\pi$ rads/sample)')
    ylabel('PSD (dB)')
    title('AR PSD Estimation, n = 500','fontsize',15)
    switch i 
        case 1
            legend('True','AR(2)')
        case 2
            legend('True','AR(4)')
        case 3
            legend('True','AR(11)')
    end
    grid on 
    grid minor
    set(gcf, 'color','w');
    hold on
end


% Question 1c
% Question 1.b
n = 10000;
rem = 500;
a = [1, -2.76, 3.81, -2.65, 0.92];
x =  filter(1,a, randn(n,1));
% discard first 500 samples to remove transient filter effects
x = x(rem+1:end);
% obtain the TRUE power spectrum 

% freqz returns th n-point freq response, as well as the angular freq for a
% digitial filter defined by the function coefficients stored in b and a
% [h,w] = freqz(b,a,n) 

[h,w] = freqz(1,a,length(x));
thePSD = abs(h).^2;
ordersOI = [2,4,6];
colors = {'r','m','c',[0.5,0,0.5],[0.5,0.5,0.5]};
MSEs1 = [];
for i = 1: length(ordersOI)
    
    % The Yule-Walker method provides an autoregressive power spectral
    % density estimate
    [aest,e] = aryule(x,ordersOI(i));
    [h,w] = freqz(1,aest,length(x));
    estPSD = abs(h).^2;
    subplot(2,3,i+3)
    plot(w/pi,10*log10(thePSD),'b','LineWidth',2);
    hold on
    plot(w/pi, 10*log10(estPSD),'r','LineWidth',2)
    MSEs1 = [MSEs1,mean((10*log10(thePSD) - 10*log10(estPSD)).^2)];
    ax = gca;
    ax.FontSize = 14;
    axis tight
    xlabel('Normalised Frequency ($\pi$ rads/sample)')
    ylabel('PSD (dB)')
    title('AR PSD Estimation, n = 9500','fontsize',15)
    switch i 
        case 1
            legend('True','AR(2)')
        case 2
            legend('True','AR(4)')
        case 3
            legend('True','AR(6)')
    end
    grid on 
    grid minor
    set(gcf, 'color','w');
    hold on
end


% 
% %% MSE
% n = 1000;
% rem = 500;
% a = [1, -2.76, 3.81, -2.65, 0.92];
% x =  filter(1,a, randn(n,1));
% % discard first 500 samples to remove transient filter effects
% x = x(rem+1:end);
% % obtain the TRUE power spectrum 
% 
% % freqz returns th n-point freq response, as well as the angular freq for a
% % digitial filter defined by the function coefficients stored in b and a
% % [h,w] = freqz(b,a,n) 
% [h,w] = freqz(1,a,length(x));
% thePSD = abs(h).^2;
% ordersOI = 1:15;
% colors = {'r','m','c',[0.5,0,0.5],[0.5,0.5,0.5]};
% MSEs0 = [];
% for i = 1: length(ordersOI)
%     
%     % The Yule-Walker method provides an autoregressive power spectral
%     % density estimate
%     [aest,e] = aryule(x,ordersOI(i));
%     [hest,w] = freqz(1,aest,length(x));
%     estPSD = abs(hest).^2;
%     MSEs0 = [MSEs0,mean((10*log10(thePSD) - 10*log10(estPSD)).^2)];
%     
% end
% 
% 
% % Question 1c
% % Question 1.b
% n = 10000;
% rem = 500;
% a = [1, -2.76, 3.81, -2.65, 0.92];
% x =  filter(1,a, randn(n,1));
% % discard first 500 samples to remove transient filter effects
% x = x(rem+1:end);
% [h,w] = freqz(1,a,length(x));
% thePSD = abs(h).^2;
% ordersOI = 1:15;
% colors = {'r','m','c',[0.5,0,0.5],[0.5,0.5,0.5]};
% MSEs1 = [];
% for i = 1: length(ordersOI)
%     
%     % The Yule-Walker method provides an autoregressive power spectral
%     % density estimate
%     [aest,e] = aryule(x,ordersOI(i));
%     [hest,w] = freqz(1,aest,length(x));
%     estPSD = abs(hest).^2;
%     MSEs1 = [MSEs1,mean((10*log10(thePSD) - 10*log10(estPSD)).^2)];
% end
