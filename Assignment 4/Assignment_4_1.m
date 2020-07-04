%% Assignment 4 1)
close all; clear all ; clc

%% Loading data
load 'time-series.mat';

%Visualising data
figure;
plot(y,'Linewidth',1)
hold on
yline(mean(y),'r--','Linewidth',1);
grid on
grid minor
title('Original time series','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]','E\{y[n]\}','Interpreter','latex')

%Removing the mean of the time-series
y_zero_mean = y - mean(y);

figure;
plot(y_zero_mean,'Linewidth',1)
hold on
yline(mean(y_zero_mean),'r--','Linewidth',1);
grid on
grid minor

%% Implementing LMS for prediction
mu=10^(-5);
num_coef = 4; %Assuming AR(4) process

[y_hat,error,weights] = LMS_pred_AR(y_zero_mean,mu,num_coef);

figure;
subplot(1,2,1)
plot(y_zero_mean,'Linewidth',1)
hold on
plot(y_hat,'Linewidth',1)
grid on
grid minor
title('One-step ahead prediction using LMS','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]-E\{y[n]\}','AR(4) estimate','Interpreter','latex')
ylim([-50 50])

subplot(1,2,2)
plot(y_zero_mean,'Linewidth',1)
hold on
plot(y_hat,'Linewidth',1)
grid on
grid minor
title('One-step ahead prediction using LMS','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]-E\{y[n]\}','AR(4) estimate','Interpreter','latex')
ylim([-50 50])
xlim([800 1000])

%% Calculating MSE and Prediction Gain
%Mean Square Error (dB)
MSE = pow2db(mean(abs(error).^2)); %In dB

figure;
plot(pow2db(abs(error).^2))

%Prediction Gain (dB)
Rp = pow2db(var(y_hat)/var(error));

%% Checking convergence
figure;
for i=1:size(weights,1)
    plot(weights(i,:))
    hold on    
end

%% Finding MSE and Rp prior and post convergence
MSE_prior = pow2db(mean(abs(error(1:200)).^2)); %In dB
Rp_prior = pow2db(var(y_hat(1:200))/var(error(1:200)));

MSE_post = pow2db(mean(abs(error(201:end)).^2)); %In dB
Rp_post= pow2db(var(y_hat(201:end))/var(error(201:end)));



