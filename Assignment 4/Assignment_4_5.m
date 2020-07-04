%% Assignment 4 5)
close all; clear all ; clc

%% Loading data
load 'time-series.mat';

%% Considering bias in prediction
a=(1:0.01:100);
mu=10^(-7);
num_coef = 4; %Assuming AR(4) process
MSE_bias=zeros(1,length(a));
Rp_bias=zeros(1,length(a));

for i=1:length(a)
    % Implementing LMS for prediction
    [y_hat,error,~] = LMS_AR_tanh_scale_bias(y,mu,num_coef,a(i));
    
    %Mean Square Error (dB)
    MSE_bias(i) = pow2db(mean(error.^2)); %In dB
    %Prediction Gain (dB)
    Rp_bias(i) = pow2db(var(y_hat)/var(error));
end

[MSE_min_bias,index_MSE_bias] = min(MSE_bias);
[Rp_max_bias,index_Rp_bias] = max(Rp_bias);

a_opt_bias = mean([a(index_MSE_bias) a(index_Rp_bias)]);

% Implementing LMS for prediction
[y_hat_bias,error_bias,weights_bias] = LMS_AR_tanh_scale_bias(y,mu,num_coef,a_opt_bias);

figure;
plot(y)
hold on
plot(y_hat_bias)

%Mean Square Error (dB)
MSE_opt_bias = pow2db(mean(abs(error_bias).^2)); %In dB
%Prediction Gain (dB)
Rp_opt_bias = pow2db(var(y_hat_bias)/var(error_bias));
    

%% Pre-training of coefficients
N_samples = 20;
N_epochs = 100;
% y_pre_training = repmat(y(1:N_samples), 1, N_epochs);
MSE_pre=zeros(1,length(a));
Rp_pre=zeros(1,length(a));
w_init_mtx = zeros(length(a),num_coef+1);
y_pre = y(1:20); %Vector used for pre-training
w_init = zeros(num_coef+1,1);
for i=1:length(a)
    %Pretraining procedure
    for epoch_ind = 1:N_epochs
        % Implementing LMS for prediction for first 20 samples
        [~,~,wpre] = LMS_AR_pretraining(y_pre,mu,num_coef,a(i),w_init); 
        w_init = wpre(:,end);
    end
    %Storing initial weights in a matrix - for optimal scale
    w_init_mtx(i,:)= w_init;
    %LMS using the initial weights found
    [y_hat_pre,error_pre,weights_pre]=LMS_AR_pretraining(y,mu,num_coef,a(i),w_init);
    
    %Mean Square Error (dB)
    MSE_pre(i) = pow2db(mean(abs(error_pre).^2)); %In dB
    %Prediction Gain (dB)
    Rp_pre(i) = pow2db(var(y_hat_pre)/var(error_pre));   
end

%Finding optimal scale
[MSE_min_pre,index_MSE_pre] = min(MSE_pre);
[Rp_max_pre,index_Rp_pre] = max(Rp_pre);

a_opt_pre = mean([a(index_MSE_pre) a(index_Rp_pre)]);

%Find corresponding initial weights
w_init_opt = w_init_mtx(round(0.5*(index_MSE_pre+index_Rp_pre)),:);

% Implementing LMS for prediction using optimal scale
[y_hat_pre_opt,error_bias,weights] = LMS_AR_pretraining(y,mu,num_coef,a_opt_pre,w_init_opt);

figure;
subplot(1,2,2)
plot(y,'Linewidth',1)
hold on
plot(y_hat_pre_opt,'Linewidth',1)
grid on
grid minor
title('LMS prediction using weight pre-training','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]','LMS with pretraining','Interpreter','latex')
ylim([-50 50])

subplot(1,2,1)
plot(y,'Linewidth',1)
hold on
plot(y_hat_pre_opt,'Linewidth',1)
grid on
grid minor
title('First 100 time-steps','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]','LMS with pretraining','Interpreter','latex')
ylim([-50 50])
xlim([0 100])

%% Observing convergence
figure;
plot(weights(1,:),'Linewidth',1)
hold on
plot(weights_bias(1,:),'Linewidth',1)
ylim([0 0.12])
grid on
grid minor
xlabel('Sample Index','Fontsize',11)
ylabel('Weight Magnitude','Fontsize',11)
title('Time-evolution of w_{0} for the two methods','Fontsize',11)
legend('With Pre-training','Without Pre-training')

%% Percentage absolute differences in mean
Perc_bias = 100*(abs(mean(y_hat_bias)-mean(y))/mean(y))
Perc_pre = 100*(abs(mean(y_hat_pre_opt)-mean(y))/mean(y))