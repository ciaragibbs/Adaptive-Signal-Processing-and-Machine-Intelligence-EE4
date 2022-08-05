%% Adaptive Signal Processing
%% 3.3 A Real Time Spectrum Analyser Using Least Mean Square
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 02/04/22
clear
close all

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

% Question 3.2a
sampNo = 1500;
nn = 1:sampNo;
fs = 2000;
var = 0.05;
eta = sqrt(var).*randn(1,sampNo) + 1j*sqrt(var).*randn(1,sampNo); % circular complex-valued white noise

% generate f(n)
phi_dot = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
% integrate dphi_dn to get phi(n)
phi = cumtrapz(phi_dot);

% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + eta;
L = 1024;
w = (0:(L-1)) .* (fs / L);
% unlike 3.2, x is now the collection of basis vectors defining the DFT
x = (1/L)*exp(1j*2*(0:sampNo-1)'*pi*(0:(L-1))/L).';
gamma = [0,0.01,0.1,0.5];
mu = 1;
count = 1;
gammas = [0,0.01,0.1,0.5];

figure
hold on
for gamma = gammas
    %pre allocation for coefficinets, error and signal
    a = complex(zeros(1,sampNo));  % need to change notation for psd code part
    e = complex(zeros(1,sampNo));

    
    %implement clms algorithm 
    [a,e] = clms_dft(x,y,mu,gamma,L);
    H = abs(a);
    %Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    subplot(2,2,count)
    surf(1:sampNo, (0:(L-1)).*(fs/L), H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time Index (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    title(sprintf('$\\gamma$ = %0.3f, $\\mu$ = 1',gamma),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    ylim([0 700])
    count = count+1;
    hold on
end

sgtitle('Time Frequency Estimation using DFT-CLMS','Interpreter','latex','fontsize',20)
set(gcf,'color','w')

%% question 3.3b
% Implementing the DFT-CLS algorithm on a segment of the EEG signal
load('EEG_Data_Assignment1.mat')
sampNo = 1200;
nn = 1:sampNo;
fs = 1200;
% generate y(n) - now it is a segment of the data!
aa = 1000;
y = POz(aa:aa+sampNo-1);
y =  y-mean(y);
L = 1024;
w = (0:(L-1)) .* (fs / L);
% unlike 3.2, x is now the collection of basis vectors defining the DFT
x = (1/L)*exp(1j*2*(0:sampNo-1)'*pi*(0:(L-1))/L).';
mu = 1;
count = 1;
gammas = [0,0.1];
figure
hold on
for gamma = gammas
    %pre allocation for coefficinets, error and signal
    a = complex(zeros(1,sampNo));  % need to change notation for psd code part
    e = complex(zeros(1,sampNo));

    
    %implement clms algorithm 
    [a,e] = clms_dft(x,y,mu,gamma,L);
    H = abs(a);
    %Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    subplot(1,2,count)
    surf(1:sampNo, (0:(L-1)).*(fs/L), H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time Index (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    title(sprintf('$\\gamma$ = %0.3f, $\\mu$ = 1',gamma),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    ylim([0 100])
    count = count+1;
    hold on
end

sgtitle('EEG Time Frequency Estimation using DFT-CLMS','Interpreter','latex','fontsize',20)
set(gcf,'color','w')
