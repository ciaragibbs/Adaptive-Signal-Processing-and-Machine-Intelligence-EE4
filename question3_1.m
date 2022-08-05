%% Adaptive Signal Processing
%% 3.1 Widely Linear Filtering and Adaptive Spectrum Estimation
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 08/04/22
clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

%% Complex LMS and Widely Linear Modelling

% Question 3.1a

sampNo = 1000;
its = 100;
% params
b = [complex(1.5,1), complex(2.5,-0.5)];
a = 1;
order = length(b);
mu = 0.01;

% pre allocation for coefficinets, error and signal
h = complex(zeros(its, order, sampNo)); 
e = complex(zeros(its, sampNo));
y = complex(zeros(1,sampNo)); % need to make the preallocation complex
ha = h;
ga = h;
ea = e;

for i=1:its
   
    x = randn(1,sampNo)+1j*randn(1,sampNo);
    % create the signal representing a first order widely linear moving
    % average process
    
    % note that the first element of the signal is 0
    y(2:sampNo) = a*x(2:sampNo) + b(1)*x(1:sampNo-1) + b(2)*conj(x(1:sampNo-1));
    
    % implement clms algorithm 
    [h(i,:,:),e(i,:)] = clms(x,y,mu,order);
    % implement aclms algorithm
    [ga(i,:,:),ha(i,:,:),ea(i,:)] = aclms(x,y,mu,order);
 

end
    
figure
hold on
subplot(1,2,1)
scatter(real(y),imag(y),20,'filled')
hold on
scatter(real(x),imag(x),20,'filled','r')
xlabel("$\Re$",'fontsize',18)
ylabel("$\Im$",'fontsize',18)
title('Circularity Evaluation','Interpreter','latex','fontsize',20)
xlim([-10 10])
ylim([-10 10])
legend('WLMA(1)','WGN')
ax = gca;
ax.FontSize = 18;
grid on
grid minor
set(gcf,'color','w')

en = abs(e).^2;
en = squeeze(mean(en, 1));
ean = abs(ea).^2;
ean = squeeze(mean(ean, 1));

subplot(1,2,2)
plot(10*log10(en),'LineWidth',2)
hold on
plot(10*log10(ean),'r','LineWidth',2)
title('Learning Curves CLMS and ACLMS Algorithms', 'fontsize', 20)
xlabel('Sample Index (n)','fontsize',18)
ylabel('Squared Error(dB)','fontsize',18)
legend('CLMS','ACLMS')
ax = gca;
ax.FontSize = 18;
grid on
grid minor

%% Question 1.b
close all
lowWindDat = load('low-wind.mat');
medWindDat = load('medium-wind.mat');
highWindDat = load('high-wind.mat');

% formulate the complex-valued wind signals
v_low = complex(lowWindDat.v_east,lowWindDat.v_north);
v_med = complex(medWindDat.v_east,medWindDat.v_north);
v_high = complex(highWindDat.v_east,highWindDat.v_north);

% determine the circularity of the respective signals
circ_low = abs(mean((v_low).^2)/mean(abs(v_low).^2));
circ_med = abs(mean((v_med).^2)/mean(abs(v_med).^2));
circ_high =  abs(mean((v_high).^2)/mean(abs(v_high).^2));

% centre of mass
com_low = [mean(real(v_low)),mean(imag(v_low))];
com_med = [mean(real(v_med)),mean(imag(v_med))];
com_high = [mean(real(v_high)),mean(imag(v_high))];

figure
% low speed
subplot(1,3,1)
scatter(real(v_low),imag(v_low),10,'filled','b')
hold on
scatter(com_low(1),com_low(2),20,'filled','y')
xlabel("$\Re$ (East)",'fontsize',16)
ylabel("$\Im$ (North)",'fontsize',16)
title(sprintf('Low-speed Data: $\\rho$ = %0.2f',circ_low),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor

% med speed
subplot(1,3,2)
scatter(real(v_med),imag(v_med),10,'filled','r')
hold on
scatter(com_med(1),com_med(2),20,'filled','y')
xlabel("$\Re$ (East)",'fontsize',16)
ylabel("$\Im$ (North)",'fontsize',16)
title(sprintf('Medium-speed Data: $\\rho$ = %0.2f',circ_med),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor

% high speed
subplot(1,3,3)
scatter(real(v_high),imag(v_high),10,'filled','m')
hold on
scatter(com_high(1),com_high(2),20,'filled','y')
xlabel("$\Re$ (East)",'fontsize',16)
ylabel("$\Im$ (North)",'fontsize',16)
title(sprintf('High-speed Data: $\\rho$ = %0.2f',circ_high),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor
set(gcf,'color','w')
%%
% Configuring the CLMS and ACLMS filters as coded in part a to fit a prediction setting

% choosing a small step size so that you account for the highest speed!
mu = [0.1;0.005;0.001];
sampNo = length(v_high); % now pre-set
filtLen = 25;
% filter length is essentially the order variable i.e. how many b terms are
% we going to have in our expression for the complex-valued signal of
% interest
all_sigs = [v_low, v_med, v_high];
figure
hold on
for j = 1:3
    
    % pre allocation error (since no iterations - no need to store h and g)
    e = complex(zeros(length(filtLen), sampNo));
    ea = e;
    sig = all_sigs(:,j);
    
    for i=1:filtLen

        % to use the clms and aclms algorithms in a 'prediction' setting
        % the observation x will be a time shifted version of the wind-speed
        % signal (back by step 1), then y will be one step ahead (the original signal basically)
        %, i.e. the prediction in the future we are aiming for
        x = [0; sig(1:end-1)]';

        % implement clms algorithm 
        [~,e(i,:)] = clms_b(x,sig',mu(j),i);
        % implement aclms algorithm
        [~,~,ea(i,:)] = aclms_b(x,sig',mu(j),i);


    end
    
    en = abs(e).^2;
    en = squeeze(mean(en, 2));
    ean = abs(ea).^2;
    ean = squeeze(mean(ean, 2));
    
    subplot(1,3,j)
    plot(10*log10(en),'LineWidth',2)
    [~,idx1] = min(10*log10(en))
    hold on
    plot(10*log10(ean),'r','LineWidth',2)
    [~,idx2] = min(10*log10(ean))
    if j == 1
    	title('Learning Curves - Low Speed', 'fontsize', 16)
    elseif j == 2
        title('Learning Curves - Mid Speed', 'fontsize', 16)
    else
        title('Learning Curves - High Speed', 'fontsize', 16)
    end
    xlabel('Filter Order')
    ylabel('MSPE (dB)')
    legend('CLMS','ACLMS')
    ax = gca;
    ax.FontSize = 16;
    grid on
    grid minor
    axis tight
    
end

set(gcf,'color','w')

%% Three Phase Power Systems: Circularity

% Question 3.1c

% a power system is balanced when:
% i) all the three signals have the same peak voltage (will just set as 1
% here)
% ii) the three voltages have a phase shift of 120 degrees, therefore
% delta_b = 0 and delta_c = 0

sampNo = 1000;
n = 1: sampNo; % discrete time vector 
fo = 100; % system frequency
fs = 10000; % sampling frequency
V = 1; % peak voltages (all match for a balanced system)
phi = [0 -2*pi/3 2*pi/3];
delta_b = 0;
delta_c = 0;

% set the voltages
v_a = V*cos(2*pi*(fo/fs)*n + phi(1));
v_b = V*cos(2*pi*(fo/fs)*n + delta_b + phi(2));
v_c = V*cos(2*pi*(fo/fs)*n + delta_c + phi(3));

% create a vector of all voltages
v_all = [v_a; v_b; v_c];

% apply the clarke transform
clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
v_projected = clarkeMat*v_all;

% need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
% j*v_beta(n)
v_clarke = complex(v_projected(2,:),v_projected(3,:));
V = [1;1;1];
% what is the circularity of the simulated system?
circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));
figure
subplot(1,2,1)
scatter(real(v_clarke), imag(v_clarke), 20,'filled');
% title(sprintf("Complex Voltage for Balanced System: $\\rho$ = %.2f" ,circularity),'interpreter','latex','fontsize',16);
title(sprintf("Complex Voltage for Balanced System: $\\rho$ = %.2f,\n $V_{a}$ = %.2f, $V_{b}$ = %.2f, $V_{c}$ = %.2f, $\\Delta_{b}$ = %.2f,$\\Delta_{c}$ = %.2f " ,circularity,V(1),V(2),V(3),delta_b,delta_c),'interpreter','latex')
xlabel("$\Re$")
ylabel("$\Im$")
ax = gca;
ax.FontSize = 16;
grid on
grid minor
set(gcf,'color','w')

% Unbalanced system - create with both alterations in magnitude and phase
% magnitudes are no longer equal = property of unbalanced system
V = [1;3;5];

% phase shifts are no longer zero = property of unbalanced system
delta_b = 0.2;
delta_c = 0.8;
phi = [0 -2*pi/3 2*pi/3];
% set the voltages
v_a = V(1)*cos(2*pi*(fo/fs)*n + phi(1));
v_b = V(2)*cos(2*pi*(fo/fs)*n + delta_b + phi(2));
v_c = V(3)*cos(2*pi*(fo/fs)*n + delta_c + phi(3));

% create a vector of all voltages
v_all = [v_a; v_b; v_c];

% apply the clarke transform
clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
v_projected = clarkeMat*v_all;

% need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
% j*v_beta(n)
v_clarke = complex(v_projected(2,:),v_projected(3,:));

% what is the circularity of the simulated system?
circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));

subplot(1,2,2)
scatter(real(v_clarke), imag(v_clarke), 20,'r','filled');
title(sprintf("Complex Voltage for Unbalanced System: $\\rho$ = %.2f,\n $V_{a}$ = %.2f, $V_{b}$ = %.2f, $V_{c}$ = %.2f, $\\Delta_{b}$ = %.2f,$\\Delta_{c}$ = %.2f " ,circularity,V(1),V(2),V(3),delta_b,delta_c),'interpreter','latex')
xlabel("$\Re$")
ylabel("$\Im$")
ax = gca;
ax.FontSize = 16;
grid on
grid minor
set(gcf,'color','w')

%% Iterate across various magnitudes

delta_b = [0;0.5;1];
delta_c = [0;1;0.5];
colors = {'b','r','m'};
V=  [1;1;1];
phi = [0 -2*pi/3 2*pi/3];
circles = [];
figure
hold on
for i = 1: length(delta_b)

    % set the voltages
    v_a = V(1)*cos(2*pi*(fo/fs)*n + phi(1));
    v_b = V(2)*cos(2*pi*(fo/fs)*n + delta_b(i) + phi(2));
    v_c = V(3)*cos(2*pi*(fo/fs)*n + delta_c(i) + phi(3));

    % create a vector of all voltages
    v_all = [v_a; v_b; v_c];

    % apply the clarke transform
    clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
    v_projected = clarkeMat*v_all;

    % need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
    % j*v_beta(n)
    v_clarke = complex(v_projected(2,:),v_projected(3,:));

    % what is the circularity of the simulated system?
    circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));
    circles = [circles, circularity];
    subplot(1,2,1)
    scatter(real(v_clarke), imag(v_clarke), 20,colors{i},'filled');
    title(sprintf("Complex Voltage for Unbalanced System: $V_{a} = V_{b} = V_{c} = 1$"), 'interpreter','latex')
    ax = gca;
    ax.FontSize = 15;
    xlabel("$\Re$");
    ylabel("$\Im$");
    grid on
    grid minor
    hold on
    
end
set(groot,'defaultLegendInterpreter','latex');
legend([strcat('$\rho$ = ',num2str(circles(1))) strcat(', $\Delta_{b}$ = ',num2str(delta_b(1))) strcat(', $\Delta_{c}$= ',num2str(delta_c(1)))],[strcat('$\rho$ = ',num2str(circles(2))) strcat(', $\Delta_{b}$ = ',num2str(delta_b(2))) strcat(', $\Delta_{c}$= ',num2str(delta_c(2)))],[strcat('$\rho$ = ',num2str(circles(3))) strcat(', $\Delta_{b}$ = ',num2str(delta_b(3))) strcat(', $\Delta_{c}$= ',num2str(delta_c(3)))] ,'Interpreter','latex')
set(gcf,'color','w')

% Iterate across various Phases

delta_b = [0;0;0];
delta_c = [0;0;0];
colors = {'b','r','m'};
V=  [1, 1,0; 1 0 1; 0 1 1];
phi = [0 -2*pi/3 2*pi/3];
circles = [];
hold on
for i = 1: length(delta_b)

    % set the voltages
    v_a = V(1,i)*cos(2*pi*(fo/fs)*n + phi(1));
    v_b = V(2,i)*cos(2*pi*(fo/fs)*n + delta_b(i) + phi(2));
    v_c = V(3,i)*cos(2*pi*(fo/fs)*n + delta_c(i) + phi(3));

    % create a vector of all voltages
    v_all = [v_a; v_b; v_c];

    % apply the clarke transform
    clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
    v_projected = clarkeMat*v_all;

    % need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
    % j*v_beta(n)
    v_clarke = complex(v_projected(2,:),v_projected(3,:));

    % what is the circularity of the simulated system?
    circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));
    circles = [circles, circularity];
    subplot(1,2,2)
    scatter(real(v_clarke), imag(v_clarke), 20,colors{i},'filled');
    title(sprintf("Complex Voltage for Unbalanced System: $\\Delta_{b} = \\Delta_{c} = 0$"), 'interpreter','latex')
    xlabel("$\Re$")
    ylabel("$\Im$")
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    hold on
    
end
legend([strcat('$\rho$ = ',num2str(circles(1))) strcat(', $V_{a}$ = ',num2str(V(1,1))) strcat(', $V_{b}$= ',num2str(V(2,1))) strcat(', $V_{c}$= ',num2str(V(3,1)))],[strcat('$\rho$ = ',num2str(circles(2))) strcat(', $V_{a}$ = ',num2str(V(1,2))) strcat(', $V_{b}$= ',num2str(V(2,2))) strcat(', $V_{c}$= ',num2str(V(3,2)))],[strcat('$\rho$ = ',num2str(circles(3))) strcat(', $V_{a}$ = ',num2str(V(1,3))) strcat(', $V_{b}$= ',num2str(V(2,3))) strcat(', $V_{c}$= ',num2str(V(3,3)))], 'Interpreter','latex')
set(gcf,'color','w')

%% Question 3.1e
close all
sampNo = 300;
n = 1: sampNo; % discrete time vector 
fo = 50; % system frequency
fs = 1000; % sampling frequency
V = 1; % peak voltages (all match for a balanced system)
phi = [0 -2*pi/3 2*pi/3];
delta_b = 0;
delta_c = 0;
mu = 0.02;
its = 1;
order = 1;
% pre allocation for coefficinets, error and signal
h = complex(zeros(its, order, sampNo)); 
e = complex(zeros(its, sampNo));
y = complex(zeros(1,sampNo)); % need to make the preallocation complex
ha = h;
ga = h;
ea = e;

% set the voltages
v_a = V*cos(2*pi*(fo/fs)*n + phi(1));
v_b = V*cos(2*pi*(fo/fs)*n + delta_b + phi(2));
v_c = V*cos(2*pi*(fo/fs)*n + delta_c + phi(3));

% create a vector of all voltages
v_all = [v_a; v_b; v_c];

% apply the clarke transform
clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
v_projected = clarkeMat*v_all;

% need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
% j*v_beta(n)
v_clarke = complex(v_projected(2,:),v_projected(3,:));

% use complex least mean squares
[h, err] = clms_e(v_clarke,mu, 1);
% what is the system frequency?
fo_pred1 = (fs/(2*pi))*atan(imag(h)./ real(h));
% use augmented complex mean squares
[g,h,ea] = aclms_e(v_clarke,mu,1);
% what is the system frequency?
fo_pred2 = ((fs/(2*pi))*atan(sqrt((imag(h)).^2 - abs(g).^2)./real(h)));

% what is the circularity of the simulated system?
circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));
figure
hold on
subplot(1,2,1)
plot(abs(fo_pred1),'b','LineWidth',2)
hold on
plot(abs(fo_pred2),'r','LineWidth',2)
hold on
yline(fo,'k--','LineWidth',2)
title(sprintf("Frequency Estimate for Balanced System: $\\rho$ = %.2f" ,circularity),'interpreter','latex','fontsize',14);
xlabel("Sample Index (n)")
ylabel('Frequency fo')
ylim([0 100])
legend('fo - CLMS','fo - ACLMS','fo - True Value','Interpreter','latex')
ax = gca;
ax.FontSize = 16;
grid on
grid minor
set(gcf,'color','w')

% repeat for unbalanced system


% magnitudes are no longer equal = property of unbalanced system
V = [1;3;5];

% phase shifts are no longer zero = property of unbalanced system
delta_b = 0.2;
delta_c = 0.8;
phi = [0 -2*pi/3 2*pi/3];
% set the voltages
v_a = V(1)*cos(2*pi*(fo/fs)*n + phi(1));
v_b = V(2)*cos(2*pi*(fo/fs)*n + delta_b + phi(2));
v_c = V(3)*cos(2*pi*(fo/fs)*n + delta_c + phi(3));

% create a vector of all voltages
v_all = [v_a; v_b; v_c];

% apply the clarke transform
clarkeMat =  sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ];
v_projected = clarkeMat*v_all;

% need to convert to a complex Clarke voltage i.e. v(n) = v_alpha(n) +
% j*v_beta(n)
v_clarke = complex(v_projected(2,:),v_projected(3,:));

% what is the circularity of the simulated system?
circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));

% use complex least mean squares
[h, err] = clms_e(v_clarke,mu, 1);
% what is the system frequency?
fo_pred1 = (fs/(2*pi))*atan(imag(h)./ real(h));
% use augmented complex mean squares
[g,h,ea] = aclms_e(v_clarke,mu,1);
% what is the system frequency?
fo_pred2 = ((fs/(2*pi))*atan(sqrt((imag(h)).^2 - abs(g).^2)./real(h)));

% what is the circularity of the simulated system?
circularity = abs(mean((v_clarke).^2)/mean(abs(v_clarke).^2));
subplot(1,2,2)
plot(abs(fo_pred1),'b','LineWidth',2)
hold on
plot(abs(fo_pred2),'r','LineWidth',2)
hold on
yline(fo,'k--','LineWidth',2)
title(sprintf("Frequency Estimate for Unbalanced System: $\\rho$ = %.2f" ,circularity),'interpreter','latex','fontsize',15);
xlabel("Sample Index (n)")
ylabel('Frequency fo')
legend('fo - CLMS','fo - ACLMS','fo - True Value','Interpreter','latex')
ylim([0 100])
ax = gca;
ax.FontSize = 16;
grid on
grid minor
set(gcf,'color','w')