%% Assignment 4 1)
close all; clear all ; clc

%% Loading data
load 'time-series.mat';

%Removing the mean of the time-series
y_zero_mean = y - mean(y);

%% Implementing LMS for prediction
mu=10^(-5);
num_coef = 4; %Assuming AR(4) process

[y_hat,error,weights] =LMS_AR_tanh(y_zero_mean,mu,num_coef);

figure;
subplot(1,2,1)
plot(y_zero_mean,'Linewidth',1)
hold on
plot(y_hat,'Linewidth',1)
grid on
grid minor
title('LMS prediction using tanh activation function','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]-E\{y[n]\}','LMS prediction with tanh','Interpreter','latex')
ylim([-50 50])

subplot(1,2,2)
plot(y_zero_mean,'Linewidth',1)
hold on
plot(y_hat,'Linewidth',1)
grid on
grid minor
title('Final 200 time-steps','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]-E\{y[n]\}','LMS prediction with tanh','Interpreter','latex')
ylim([-50 50])
xlim([800 1000])

%% Calculating MSE and Prediction Gain
%Mean Square Error (dB)
MSE = pow2db(mean(abs(error).^2)); %In dB

figure;
plot(pow2db(abs(error).^2))

%Prediction Gain (dB)
Rp = pow2db(var(y_hat)/var(error));

