%% Q4 From LMS to Deep Learning
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 03/04/22
clear
close all
clc
load('time-series.mat')

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

%% Question 4.1
load('time-series.mat')
sampNo = length(y);
n = 1: sampNo;
mu = 0.00001;
gamma = 0;
order = 4;
y = y -mean(y);

[yhat,w,error] = LMS(y,mu,gamma,order);

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',25)
legend('Centred Time Series','Estimated Time Series')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(750:1000,y(750:end),'b','LineWidth',1.5)
hold on
plot(750:1000,yhat(750:end),'r','LineWidth',1.5)
xlim([750,1000])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',25)
legend('Centred Time Series','Estimated Time Series')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));


MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE_end);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));


%% Question 4.2

% building upon the LMS algorithm to add a dynamical perceptron
sampNo = length(y);
n = 1: sampNo;
mu = 0.00001;
gamma = 0;
order = 4;
y = y -mean(y);
alpha = 1;

[yhat,w,error] = LMS_dypn(y,mu,order,alpha,0,[0;0;0;0],0);

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Zero-Mean Time Series','Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(750:1000,y(750:end),'b','LineWidth',1.5)
hold on
plot(750:1000,yhat(750:end),'r','LineWidth',1.5)
xlim([750,1000])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Zero-Mean Time Series','Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));


MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE_end);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));

%% Question 4.3
load('time-series.mat')
% Selecting an optimal value of alpha, ranging between 40 and 100 for our guess
% building upon the LMS algorithm to add a dynamical perceptron
sampNo = length(y);
n = 1: sampNo;
mu = 0.0000001;
gamma = 0;
order = 4;
y = y -mean(y);
%y = [ones(1,order),y']';
starter =  40;
step = 0.1;
ender = 100;
alphas = [starter:step:ender];
MSEs= [];
R_ps = [];
for alpha = alphas
    
    [yhat,w,error] = LMS_dypn(y,mu,order,alpha,0,[0;0;0;0],0);
    % calculating the MSE between the true and estimated signals
    MSEs = [MSEs, 10*log10(mean(abs(error).^2))];
    %MSE_db = 10*log10(MSE);
    R_ps = [R_ps,10*log10(var(yhat)/var(error))];
end

figure
subplot(1,2,1)
plot(alphas,MSEs,'b','LineWidth',2)
hold on
[val,ind] = min(MSEs);
plot(starter + ind*step,val,'r*','MarkerSize',10)
xlabel('$\alpha$','fontsize',18)
ylabel('Mean Squared Error (dB)','fontsize',18)
ax = gca;
ax.FontSize = 15; 
grid on 
grid minor
subplot(1,2,2)
plot(alphas,R_ps,'b','LineWidth',2)
hold on
[val,ind] = max(R_ps);
plot(starter + ind*step,val,'r*','MarkerSize',10)
xlabel('$\alpha$','fontsize',18)
ylabel('Prediction Gain(dB)','fontsize',18)
ax = gca;
ax.FontSize = 15; 
grid on 
grid minor
sgtitle('Finding the Optimal Value for $\alpha$','Interpreter','Latex','fontsize',20)
set(gcf,'color','w')

% calculating the MSE between the true and estimated signals
MSEs = [MSEs, 10*log10(mean(abs(error).^2))];
%MSE_db = 10*log10(MSE);
R_ps = [R_ps,10*log10(var(yhat)/var(error))];



MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE_end);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Zero-Mean Time Series','Scaled Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(750:1000,y(750:1000),'b','LineWidth',1.5)
hold on
plot(750:1000,yhat(750:1000),'r','LineWidth',1.5)
xlim([750,1000])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',16)
legend('True Zero-Mean Time Series','Scaled Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')
%% NOTE THIS SECTION WAS NOT SHOWN IN THE COURSEWORK ONLY BRIEFLY DISCUSSED
% using the gradient update!
% Selecting an optimal value of alpha, ranging between 40 and 100 for our guess
% building upon the LMS algorithm to add a dynamical perceptron
sampNo = length(y);
n = 1: sampNo;
mu = 0.0000001;
gamma = 0;
order = 4;
y = y -mean(y);
alpha = 92.5;
MSEs= [];
R_ps = [];

[yhat,w,error,alphas] = LMS_dypn(y,mu,order,alpha,0,[0;0;0;0]);
% calculating the MSE between the true and estimated signals
MSEs = [MSEs, 10*log10(mean(abs(error).^2))];
%MSE_db = 10*log10(MSE);
R_ps = [R_ps,10*log10(var(yhat)/var(error))];

MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE_end);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Zero-Mean Time Series','Scaled Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(750:1000,y(750:end),'b','LineWidth',1.5)
hold on
plot(750:1000,yhat(750:end),'r','LineWidth',1.5)
xlim([750,1000])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Zero-Mean Time Series','Scaled Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

%% Question 4.4

% Considering the original time-series now with the non-zero mean i.e. a
% bias
% adding bias to the input model

load('time-series.mat')
% Selecting an optimal value of alpha, ranging between 40 and 100 for our guess
% building upon the LMS algorithm to add a dynamical perceptron
sampNo = length(y);
n = 1: sampNo;
mu = 0.0000001;
gamma = 0;
order = 4;
starter =  40;
step = 0.1;
ender = 100;
alpha = 75.7;
MSEs= [];
R_ps = [];
bias = 1;
winit = zeros(order+bias,1);

[yhat,w,error] = LMS_dypn(y,mu,order,alpha,0,winit,bias);

% calculating the MSE between the true and estimated signals
MSE = 10*log10(mean(abs(error).^2));
MSE_db = 10*log10(MSE);
R_ps = [R_ps,10*log10(var(yhat)/var(error))];

MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE_end);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Time Series','Scaled Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(750:1000,y(750:end),'b','LineWidth',1.5)
hold on
plot(750:1000,yhat(750:end),'r','LineWidth',1.5)
xlim([750,1000])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Time Series','Scaled Dynamical Perceptron Estimated Time Series','Location','southeast')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')
% plot the evolution of the weights
figure
plot(w(1,:),'b','LineWidth',1.5)
hold on
plot(w(2,:),'r','LineWidth',1.5)
hold on
plot(w(3,:),'m','LineWidth',1.5)
hold on
plot(w(4,:),'c','LineWidth',1.5)
hold on
plot(w(5,:),'Color',[0.5,0,0.5],'LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('Evolution of LMS Algorithm Weights (Scaled, Biased Dynamical Perceptron)','Interpreter','latex','fontsize',18)
legend('w0','w1','w2','w3','w4','w5','Location','southeast')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

%% Question 4.5

% pre-train the weights by over-fitting to a small number of samples
% starting with w(0) = 0 and using 100 iterations to fir the first 20
% samples t yield w_init

% then use w_init to predict the entire time series

load('time-series.mat')
% Selecting an optimal value of alpha, ranging between 40 and 100 for our guess
% building upon the LMS algorithm to add a dynamical perceptron
sampNo = length(y);
n = 1: sampNo;
mu = 0.0000001;
gamma = 0;
order = 4;
starter =  40;
step = 0.1;
ender = 100;
alphas = [starter:step:ender];
MSEs= [];
R_ps = [];
step = mu; % overwritting for the pretraining part..!
b = 1;
segLength = 20;
epochs = 100;
winit = zeros(order+b,1); % the first initialisation will just be winit = [0 0 0 0]'
all_winits = zeros(order+b,length(alphas));
count = 1;
MSEs = [];
R_s = [];
for alpha = alphas
    
    % training on the first 20 samples with the added bias
    for e = 1: epochs
        y_sample = y(1:segLength); % accounting for the ones added for bias
        
        % to pre-train
            xOut = zeros(length(y_sample),1);
            err = xOut;
            w = zeros(order+b,length(y_sample)+1); % since order can be > 1 (here 2)
            w(:,1) = winit;
            xShift = zeros(order,length(y_sample)); % x(n-k)
            % creating two shifted vectors by i, i.e. the order length
            for i = 1: order
                xShift(i,:) = [ zeros(1,i), y_sample(1: length(y_sample)-i)']; 
            end
            if b
                xShift = [ones(1,length(y_sample)); xShift];
            end

            for k = 1: length(y_sample)

                % calculate the prediction
                xOut(k) = alpha*tanh(w(:,k)'*xShift(:,k));
                err(k) = y_sample(k)-xOut(k);
                act_function = alpha*(1-(xOut(k)/alpha)^2);
                % update
                w(:,k+1)=w(:,k)+(step*act_function*err(k)).*xShift(:,k);
  
            end
            w =  w(:,2:end);
            winit = w(:,end);
    end
    
    all_winits(:,count) = winit;
    
    % now going ahead with the normal algorithm
    [yhat,w,error] = LMS_dypn(y,mu,order,alpha,0,winit,b);
    % calculating the MSE between the true and estimated signals
    MSEs = [MSEs, 10*log10(mean(abs(error(order+1:end)).^2))];
    %MSE_db = 10*log10(MSE);
    R_ps = [R_ps,10*log10(var(yhat(order+1:end))/var(error(order+1:end)))];
    
    count = count+1;
end
%%
figure
subplot(1,2,1)
step = 0.1;
plot(alphas,MSEs,'b','LineWidth',2)
hold on
[val,ind] = min(MSEs);
plot(starter + ind*step,val,'r*','MarkerSize',10)
xlabel('$\alpha$','fontsize',14)
ylabel('Mean Squared Error (dB)','fontsize',14)
grid on 
grid minor
subplot(1,2,2)
plot(alphas,R_ps,'b','LineWidth',2)
hold on
[val,ind] = max(R_ps);
plot(starter + ind*step,val,'r*','MarkerSize',10)
xlabel('$\alpha$','fontsize',14)
ylabel('Prediction Gain(dB)','fontsize',14)
grid on 
grid minor
sgtitle('Finding the Optimal Value for $\alpha$','Interpreter','Latex','fontsize',20)
set(gcf,'color','w')




%%

load('time-series.mat')
% now taking the winits @ alpha = 67.9
sampNo = length(y);
n = 1: sampNo;
mu = 0.0000001;
gamma = 0;
order = 4;
starter =  40;
step = 0.1;
ender = 100;
alpha = 67.9;
MSEs= [];
R_ps = [];


winit = all_winits(:,find(alphas == 67.9));
[yhat,w,error] = LMS_dypn(y,mu,order,alpha,0,winit,b);

% calculating the MSE between the true and estimated signals
MSE = 10*log10(mean(abs(error).^2));
%MSE_db = 10*log10(MSE);
R_ps = [R_ps,10*log10(var(yhat)/var(error))];

MSE_end = mean(error(750:1000).^2);
MSE_db_end = 10*log10(MSE);
R_p_end = 10*log10(var(yhat(750:1000))/var(error(750:1000)));

figure
plot(y,'b','LineWidth',1.5)
hold on
plot(yhat,'r','LineWidth',1.5)
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('True Centred Time Series','Dynamical Perceptron Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')

% showing that when we approcahing the end of the signal, the prediction
% improves
figure
plot(1:50,y(1:50),'b','LineWidth',1.5)
hold on
plot(1:50,yhat(1:50),'r','LineWidth',1.5)
xlim([1,50])
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('One-Step-Ahead Prediction of AR(4) Time Series','Interpreter','latex','fontsize',18)
legend('Centred Time Series','Estimated Time Series')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor
set(gcf,'color','w')