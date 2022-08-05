%% Adaptive Signal Processing
%% 3.2 Adaptive AR Model Based Time-Frequency Estimation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 01/04/22
clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

% Question 3.2a
sampNo = 1500;
n = 1:sampNo;
fs = 2000;
var = 0.05;
eta = sqrt(var).*randn(1,sampNo) + 1j*sqrt(var).*randn(1,sampNo); % circular complex-valued white noise

% generate f(n)
phi_dot = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
% integrate dphi_dn to get phi(n)
phi = cumtrapz(phi_dot);

% plot the f(n) and phi(n)
figure
subplot(1,2,1)
plot(phi_dot,'r','LineWidth',2)
xlabel('Time Index (n)','fontsize',18)
ylabel('Frequency (Hz)','fontsize',18)
title('Frequency $f(n)$','Interpreter','latex','fontsize',20)
grid on
grid minor

subplot(1,2,2)
plot(wrapTo2Pi(phi),'r','LineWidth',2)
xlabel('Time Index (n)','fontsize',18)
xlim([0 200])
ylabel('Angle (rads)','fontsize',18)
title('Phase $\Phi(n)$','Interpreter','latex','fontsize',20)
ax = gca;
ax.FontSize = 15;
grid on
grid minor
set(gcf,'color','w')

%
% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + eta;

% using aryule to find the AR(1) coefficinet for the complete signal

orders = [1,5,10];
colors = {[0 0 1],[1 0 0],[1 0 1],[0.5 0 0.5]};
count = 1;
figure
hold on
for i = orders

    % using aryule to find the AR(1) coefficinet for the complete signal
    a = aryule(y, i);
    % obtain the N-point frequency response vector, and the corresponding
    % angular frequency vector w - for the digitial filter constructed
    [h,w] = freqz(1,a,sampNo,fs);
    % compute the power spectral density
   
    psd = 10*log10(abs(h).^2);
    subplot(1,3,count)
    plot(w,psd,'color',colors{count},'LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)','fontsize',10)
    title(sprintf('AR(%0.0f)',i),'Interpreter','latex','fontsize',18)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    count = count+1;
end
sgtitle('Power Spectra of Frequency Modulated Signal','fontsize',18)
set(gcf,'color','w')
%
% does the method capture frequency changes in f(n)?
% splitting into 500 segment lengths in order to plot the frequency
% estimates of each phase separately
figure
hold on
for i = 1:3
    
    % using aryule to find the AR(1) coefficinet for the complete signal
    a = aryule(y(500*(i-1)+1:500*i), 1);
    % obtain the N-point frequency response vector, and the corresponding
    % angular frequency vector w - for the digitial filter constructed
    [h,w] = freqz(1,a,sampNo/3,fs);
    % compute the power spectral density
   
    psd = 10*log10(abs(h).^2);
    subplot(1,3,i)
    plot(w,psd,'b','LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)','fontsize',10)
    title(sprintf('AR(1) Segment %0.0f',i),'Interpreter','latex','fontsize',18)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    
end
sgtitle('Power Spectra of Frequency Modulated Signal','fontsize',18)
set(gcf,'color','w')
%%
allFreqs = [];
interval = 3;
fs = 10000;
% trying to improve the estimate with smaller segments
for i = 1:1500/interval
    
    for j = 1:interval
        % using aryule to find the AR(1) coefficinet for the complete signal
        try
            a = aryule(y(interval*(i-1)+j:interval*i+(j-1)), 1);
            % obtain the N-point frequency response vector, and the corresponding
            % angular frequency vector w - for the digitial filter constructed
            [h,w] = freqz(1,a,1500/interval,fs);
            % compute the power spectral density

            psd = 10*log10(abs(h).^2);
            [~,maxInd] = max(abs(psd));
            allFreqs = [allFreqs, maxInd];
        catch
        end
    end

end
figure
plot(allFreqs,'b','LineWidth',1)
hold on
plot(phi_dot,'k--','LineWidth',2)
xlabel('Time Index (n)')
ylabel('Frequency(Hz)')
title('Windowing Data to Improve AR(1) Estimate','Interpreter','latex','fontsize',18)
legend('True Frequency Values','Windowed AR(1) Estimate')
grid on
grid minor
set(gcf,'color','w')

%% Question 3.2b
sampNo = 1500;
n = 1:sampNo;
fs = 2000;
var = 0.05;
eta = sqrt(var).*randn(1,sampNo) + 1j*sqrt(var).*randn(1,sampNo); % circular complex-valued white noise

% generate f(n)
phi_dot = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];

% integrate dphi_dn to get phi(n)
phi = cumtrapz(phi_dot);

% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + eta;
x = [0, y(1:end-1)];
mu = [0.005,0.05,0.1,0.5];
count = 1;
figure
hold on
for step = mu
    % pre allocation for coefficinets, error and signal
    a = complex(zeros(1,sampNo));  % need to change notation for psd code part
    e = complex(zeros(1,sampNo));

    % implement clms algorithm 
    [a,e] = clms_b(x,y,step,1);

    L = 1024;
    H = zeros(L, sampNo);

    % code provided in the coursework documentation
    for n=1:sampNo
        % Run complex-valued LMS algorthm to esitmate AR coefficient a_hat1(n)
        [h, w] = freqz(1, [1; -conj(a(n))], L);
        H(:, n) = abs(h).^2;
    end

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    subplot(2,2,count)
    surf(1:sampNo,((w*fs)/(2*pi)), H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time Index (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    title(sprintf('$\\mu$ = %0.3f',step),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    count = count+1;
    hold on
end

sgtitle('Time Frequency Estimation using CLMS','Interpreter','latex','fontsize',20)
set(gcf,'color','w')