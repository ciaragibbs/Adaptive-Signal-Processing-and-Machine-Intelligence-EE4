%% Robust Regression
% Author: Ciara Gibbs
% CID: 01498482
% Last edit: 12/04/22
clear 
close all 
clc

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

%% Q 1.6a

load("PCAPCR.mat")

% Obtain the singular value decomposition of X and Xnoise

% SVD of X
[ux,sx,vx] = svd(X);
[un,sn,vn] = svd(Xnoise);
figure
subplot(1,2,1)
stem(sum(sx),'color','red','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel("Subspace Dimension Index ",'fontsize',15)
ylabel("Singular Value",'fontsize',15)
title("Singular Values for SVD of $X$",'Interpreter','latex','fontsize',15)
grid on
grid minor
subplot(1,2,2)
stem(sum(sn),'color','red','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel("Subspace Dimension Index ",'fontsize',15)
ylabel("Singular Value",'fontsize',15)
title("Singular Values for SVD of $X_{noise}$",'Interpreter','latex','fontsize',15)
grid on 
grid minor
set(gcf,'color','w')
%%
% Plot the square error between each singular value of X and Xnoise.
% Explain the effect on the singular values, and state at what point would
% it become hard to identify the rank of the matrix

se = sum((sx-sn).^2);
figure
subplot(1,2,1)
stem(se,'color','red','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel('Subspace Dimension Index','fontsize',15)
ylabel('Squared Error','fontsize',15)
title('Squared Error between Singular Values of $X$ and $X_{noise}$','Interpreter','latex','fontsize',14)
grid on
grid minor
set(gcf,'color','w')

%% Question 1.6b

% Using only the r most significant principal components (as determined by
% the identified rank), create a low rank approximation of Xnoise denoted
% by X~

%s = svds(A,k) returns the k largest singular values.
[usub,ssub,vsub] = svds(Xnoise,3); % 3 since this is needed for low rank approx
Xnoise_reconstruct = usub*ssub*vsub';
se1 = sum((X-Xnoise).^2);
se2 = sum((X-Xnoise_reconstruct).^2);

subplot(1,2,2)
stem(se1,'color','blue','LineWidth',2)
hold on
stem(se2,'color','red','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel('Subspace Dimension Index','fontsize',15)
ylabel('Squared Error','fontsize',15)
title('Squared Error between X, $X_{noise}$ and low rank approximation of $X_{noise}$','Interpreter','latex','fontsize',14)
legend('$X-X_{noise}$','$X-X_{reconstructed-noise}$','Interpreter','latex')
grid on 
grid minor
set(gcf,'color','w')


%% Question 1.6c

% Ordinary Least Squares Method for Matrix B

% Since X is sub-rank, the OLS solution is not achievable.
% However, Xnoise is full-rank, permiting an OLS solution - but may
% introduce correlation

B_ols = inv(Xnoise'*Xnoise)*Xnoise'*Y;
Y_ols = Xnoise*B_ols;

% Principal Component Regression Method for Matrix B
% First applies PCA
% Avoids the issue of collinearity and noise in the input matrix
B_pcr = vsub*inv(ssub)*usub'*Y;
Y_pcr = Xnoise_reconstruct*B_pcr;


error_ols =  sum((Y-Y_ols).^2);
error_pcr =  sum((Y-Y_pcr).^2);

error_ols_test = sum((Ytest-Y_ols).^2);
error_pcr_test= sum((Ytest-Y_pcr).^2);

figure

subplot(1,3,1)
stem(error_ols,'color','red','LineWidth',2)
hold on
stem(error_pcr,'color','blue','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel("Subspace Dimension Index ",'fontsize',15)
ylabel("Squared Error",'fontsize',15)
title("In-sample Squared Error",'fontsize',15)
legend('OLS','PCR','Interpreter','latex')
grid on
grid minor
subplot(1,3,2)
stem(error_ols_test,'color','red','LineWidth',2)
hold on
stem(error_pcr_test,'color','blue','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel("Subspace Dimension Index ",'fontsize',15)
ylabel("Squared Error",'fontsize',15)
title("Out-sample Squared Error",'Interpreter','latex','fontsize',15)
legend('OLS','PCR','Interpreter','latex')
grid on
grid minor

% Question 1.6d

% Best way to assess effectiveness of PCR vs OLS is via testing the
% regression coefficients B_hat over a test data set.
% regval - new realisation of test data + estimate Y_hat
% Input is the regression coefficients


err_ols = [];
err_pcr = [];

for i = 1:10
    
    [Y_ols, Y_olsh] = regval(B_ols);
    [Y_pcr, Y_pcrh] = regval(B_pcr);
    
    err_ols= [err_ols; sum((Y_ols - Y_olsh).^2)];
    err_pcr = [err_pcr; sum((Y_pcr - Y_pcrh).^2)];
    
end

subplot(1,3,3)
stem(mean(err_ols),'color','red','LineWidth',2)
hold on
stem(mean(err_pcr),'color','blue','LineWidth',2)
axis tight
ax = gca;
ax.FontSize = 15;
xlabel('Subspace Dimension Index','fontsize',15)
ylabel('Mean Squared Error','fontsize',15)
title('Mean Squared Error of Data Ensembles ','fontsize',15)
legend('OLS','PCR','Interpreter','latex')
grid on 
grid minor
set(gcf,'color','w')

