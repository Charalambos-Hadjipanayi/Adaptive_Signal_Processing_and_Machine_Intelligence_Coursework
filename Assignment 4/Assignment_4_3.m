%% Assignment 4 3)
close all; clear all ; clc

%% Loading data
load 'time-series.mat';

%Removing the mean of the time-series
y_zero_mean = y - mean(y);

%% Calculating MSE and Prediction Gain for different a
a=(1:0.05:100);
MSE=zeros(1,length(a));
Rp=zeros(1,length(a));

mu=10^(-7);
num_coef = 4; %Assuming AR(4) process

for i=1:length(a)
    
    % Implementing LMS for prediction
    [y_hat,error,~] =LMS_AR_tanh_scale(y_zero_mean,mu,num_coef,a(i));
    
    %Mean Square Error (dB)
    MSE(i) = pow2db(mean(abs(error).^2)); %In dB
    %Prediction Gain (dB)
    Rp(i) = pow2db(var(y_hat)/var(error));
end

%Finding optimal values of a
[MSE_min,index_MSE] = min(MSE);
[Rp_max,index_Rp] = max(Rp);

figure;
subplot(1,2,1)
plot(a,MSE,'Linewidth',1)
hold on
stem(a(index_MSE),MSE_min,'Linewidth',1,'Linestyle','none')
xlabel('Scale constant a','Fontsize',11)
ylabel('MSE (dB)','Fontsize',11)
title('Non-constant overall learning rate','Fontsize',11)
grid on
grid minor
xlim([1 a(end)])
ylim([5 24])

subplot(1,2,2)
plot(a,Rp,'Linewidth',1)
hold on
stem(a(index_Rp),Rp_max,'Linewidth',1,'Linestyle','none')
xlabel('Scale constant a','Fontsize',11)
ylabel('R_{p} (dB)','Fontsize',11)
title('Non-constant overall learning rate','Fontsize',11)
grid on
grid minor
xlim([1 a(end)])

a_opt_vec(1) = a(index_MSE);
a_opt_vec(2) = a(index_Rp);

a_opt = mean(a_opt_vec);

% Implementing LMS for prediction
[y_hat,error,weights] = LMS_AR_tanh_scale(y_zero_mean,mu,num_coef,a_opt);

%Mean Square Error (dB)
MSE_opt = pow2db(mean(abs(error).^2)); %In dB
%Prediction Gain (dB)
Rp_opt = pow2db(var(y_hat)/var(error));

figure;
subplot(1,2,1)
plot(y_zero_mean,'Linewidth',1)
hold on
plot(y_hat,'Linewidth',1)
grid on
grid minor
title('LMS prediction using scaled tanh activation function','Fontsize',11)
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
title('Final 200 time-steps','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]-E\{y[n]\}','LMS prediction with a\dot tanh','Interpreter','latex')
ylim([-50 50])
xlim([800 1000])

%% Observing convergence
figure;
for i=1:size(weights,1)
    plot(weights(i,:))
    hold on
end