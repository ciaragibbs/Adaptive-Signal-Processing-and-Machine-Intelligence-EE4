%% Adaptive Signal Processing
%% 2.2 The Least Mean Squares Algorithm
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 11/04/22

clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');


%% Adaptive Step Sizes

% Question 2.2a
sampNo =  1000;
its = 1000;
mu_init = 0.1;
rho = 0.001;
mus = [0.01,0.1,0.2,0.2,0.2];
alpha = 0.8;
% define the MA model
b = [1,0.9];
order = length(b) -1;
a = 1;
v = 0.5;
testNo = 5;
eta = zeros(its,sampNo);
for i = 1:its
eta(i,:) = sqrt(v)*randn(1,sampNo);
end
figure
hold on
ax = gca;
ax.FontSize = 15;
colors = {'b','r','m','c',[0.5,0,0.5]};
subplot(1,3,1)
for i = 1:5
    
    w_all = zeros(its,sampNo);
    w_error = w_all;
    xhat_all = w_all;
    MSPE = zeros(its,1);
    
    for j = 1: its
       
       x = filter(b,a,eta(j,:));
       [xHat,w,err,endMus] = LMS_GASS(eta(j,:),x,mus(i),rho,alpha,order,i);
       w_all(j,:) = w; 
    end
      w_error = b(2)*ones(its,sampNo) - w_all;
      plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
      hold on
end
ax = gca;
ax.FontSize = 15;
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',12);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0.2$','interpreter','latex')
grid on
grid minor

mus = [0.01,0.1,0.1,0.1,0.1];
subplot(1,3,2)
for i = 1:5
    
    w_all = zeros(its,sampNo);
    w_error = w_all;
    xhat_all = w_all;
    MSPE = zeros(its,1);
    
    for j = 1: its
       
       x = filter(b,a,eta(j,:));
       [xHat,w,err,endMus] = LMS_GASS(eta(j,:),x,mus(i),rho,alpha,order,i);
       w_all(j,:) = w; 
    end
      w_error = b(2)*ones(its,sampNo) - w_all;
      plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
      hold on
end
ax = gca;
ax.FontSize = 15;
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',13);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0.1$','interpreter','latex')
grid on
grid minor

mus = [0.01,0.1,0,0,0];
subplot(1,3,3)
for i = 1:5
    
    w_all = zeros(its,sampNo);
    w_error = w_all;
    xhat_all = w_all;
    MSPE = zeros(its,1);
    
    for j = 1: its
       
       x = filter(b,a,eta(j,:));
       [xHat,w,err,endMus] = LMS_GASS(eta(j,:),x,mus(i),rho,alpha,order,i);
       w_all(j,:) = w; 
    end
      w_error = b(2)*ones(its,sampNo) - w_all;
      plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
      hold on
end
ax = gca;
ax.FontSize = 15;
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',12);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0$','interpreter','latex')
grid on
grid minor
sgtitle('Weight Error Curves','fontsize',16)
set(gcf,'color','w')

%%
mus = [0.01,0.1,0.1,0.1,0.1];
figure
hold on 
for i = 1:5
    
    w_all = zeros(its,sampNo);
    w_error = w_all;
    xhat_all = w_all;
    MSPE = zeros(its,1);
    
    for j = 1: its
       
       x = filter(b,a,eta(j,:));
       [xHat,w,err,endMus] = LMS_GASS(eta(j,:),x,mus(i),rho,alpha,order,i);
       w_all(j,:) = w; 
    end
      w_error = b(2)*ones(its,sampNo) - w_all;
     
      plot(10*log10(mean(w_error.^2)),'color',colors{i},'LineWidth',1.5);
      hold on
end
ax = gca;
ax.FontSize = 15;
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',12);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Squared Weight Error (dB)','fontsize',15)
title('Squared Weight Error Curves - $\mu_{initGASS} = 0.1$','interpreter','latex')
grid on
grid minor
set(gcf,'color','w')

%% Question 2.3c

sampNo =  1000;
its = 1000;
mu_init = 0.1;
rho = 0.001;
alpha = 0.8;
% define the MA model
b = [1,0.9];
order = length(b) -1;
a = 1;
v = 0.5;
testNo = 5;
eta = zeros(its,sampNo);
for i = 1:its
eta(i,:) = sqrt(v)*randn(1,sampNo);
end
sampNo = sampNo/2;
w_allB = zeros(its,sampNo);
w_errB = w_allB;
w_allG = w_allB;
w_errG = w_allB;
mu = 0.1;
epsilon_init = 1/mu;

for j = 1: its

   x = filter(b,a,eta(j,:));
   x =  x(501:end);
   [~,wB,errB,~] = LMS_GASS(eta(j,501:end),x,0.1,0.002,alpha,order,i);
   [~,wG,errG] = GNGD(eta(j,501:end),x,1,epsilon_init,order,0.005);
   w_errB(j,:) = errB;
   w_allB(j,:) = wB; 
   w_allG(j,:) = wG;
   w_errG(j,:) = errG;
end

figure
subplot(1,3,1)
hold on
w_errorB = b(2)*ones(its,sampNo) - w_allB;
w_errorG = b(2)*ones(its,sampNo) - w_allG;
plot(nanmean(w_errorB),'b','LineWidth',1.5);
hold on
plot(nanmean(w_errorG),'r','LineWidth',1.5);
ax = gca;
ax.FontSize = 15;
legend('Benveniste','GNGD','Interpreter','latex','fontsize',13);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('Weight Error Curves','interpreter','latex')
grid on
grid minor
subplot(1,3,2)
plot(nanmean(w_errorB),'b','LineWidth',1.5);
hold on
plot(nanmean(w_errorG),'r','LineWidth',1.5);
ax = gca;
ax.FontSize = 15;
legend('Benveniste','GNGD','Interpreter','latex','fontsize',13);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
xlim([0 100])
title('Weight Error Curves (Zoom In)','interpreter','latex')
grid on
grid minor
subplot(1,3,3)
plot(nanmean(10*log10(w_errorB.^2)),'b','LineWidth',1.5);
hold on
plot(nanmean(10*log10(w_errorG.^2)),'r','LineWidth',1.5);
ax = gca;
ax.FontSize = 15;
legend('Benveniste','GNGD','Interpreter','latex','fontsize',13);
xlabel('Sample Index (n)','fontsize',15)
ylabel('Squared Weight Error (dB)','fontsize',15)
xlim([0 100])
title('Squared Weight Error Curves','interpreter','latex')
grid on
grid minor
set(gcf, 'color','w')
% %% optimisation
% sampNo =  1000;
% its = 1000;
% mu_init = 0.1;
% rho = 0.001;
% alpha = 0.8;
% % define the MA model
% b = [1,0.9];
% order = length(b) -1;
% a = 1;
% v = 0.5;
% testNo = 5;
% eta = zeros(its,sampNo);
% for i = 1:its
% eta(i,:) = sqrt(v)*randn(1,sampNo);
% end
% 
% w_allB = zeros(its,sampNo);
% w_errorB = w_allB;
% w_allG = w_allB;
% w_errorG = w_allB;
% mu = 0.1;
% epsilon_init = 1/mu;
% figure
% hold on
% for i = [0.1,0.2]
%     w_allB = zeros(its,sampNo);
%     w_errorB = w_allB;
%     for j = 1: its
% 
%        x = filter(b,a,eta(j,:));
%        [~,wB,~,~] = LMS_GASS(eta(j,:),x,i,0.001*(0.2/i),alpha,order,3);
%        %[~,wG,~] = GNGD(eta(j,:),x,1,epsilon_init,order,0.005);
%        w_allB(j,:) = wB; 
%        
%        %w_allG(j,:) = wG;
%     end
%     w_errorB = b(2)*ones(its,sampNo) - w_allB;
%     w_errorB(any(isnan(w_errorB), 1), :) = [];
%     plot(mean(w_errorB),'LineWidth',1.5);
%     hold on
% end
