%% Assignment 4 4)
close all; clear all ; clc

%% Loading data
load 'time-series.mat';

%% Calculating MSE and Prediction Gain for different a
a=(1:0.01:100);
MSE=zeros(1,length(a));
Rp=zeros(1,length(a));

mu=10^(-7);
num_coef = 4; %Assuming AR(4) process

for i=1:length(a)
    
    % Implementing LMS for prediction
    [y_hat,error,~] =LMS_AR_tanh_scale(y,mu,num_coef,a(i));
    
    %Mean Square Error (dB)
    MSE(i) = pow2db(mean(abs(error).^2)); %In dB
    %Prediction Gain (dB)
    Rp(i) = pow2db(var(y_hat)/var(error));
end

figure;
subplot(1,2,1)
plot(a,MSE)
subplot(1,2,2)
plot(a,Rp)

[MSE_min,index_MSE] = min(MSE);
[Rp_max,index_Rp] = max(Rp);

a_opt = mean([a(index_MSE) a(index_Rp)]);

% Implementing LMS for prediction
[y_hat_opt,error_opt,weights_opt] =LMS_AR_tanh_scale(y,mu,num_coef,a_opt);

figure;
plot(y)
hold on
plot(y_hat_opt)

%Mean Square Error (dB)
MSE_opt = pow2db(mean(abs(error_opt).^2)); %In dB
%Prediction Gain (dB)
Rp_opt = pow2db(var(y_hat_opt)/var(error_opt));

% Finding MSE and Rp prior and post convergence
MSE_prior_opt = pow2db(mean(abs(error_opt(1:200)).^2)); %In dB
Rp_prior_opt = pow2db(var(y_hat_opt(1:200))/var(error_opt(1:200)));

MSE_post_opt = pow2db(mean(abs(error_opt(201:end)).^2)); %In dB
Rp_post_opt = pow2db(var(y_hat_opt(201:end))/var(error_opt(201:end)));
    
%% Considering bias in prediction
MSE_bias=zeros(1,length(a));
Rp_bias=zeros(1,length(a));

for i=1:length(a)
    
    % Implementing LMS for prediction
    [y_hat,error,w] =LMS_AR_tanh_scale_bias(y,mu,num_coef,a(i));
    
    %Mean Square Error (dB)
    MSE_bias(i) = pow2db(mean(abs(error).^2)); %In dB
    %Prediction Gain (dB)
    Rp_bias(i) = pow2db(var(y_hat)/var(error));
end

[MSE_min_bias,index_MSE_bias] = min(MSE_bias);
[Rp_max_bias,index_Rp_bias] = max(Rp_bias);

a_opt_bias = mean([a(index_MSE_bias) a(index_Rp_bias)]);

% Implementing LMS for prediction
[y_hat_bias,error_bias,weights_bias] = LMS_AR_tanh_scale_bias(y,mu,num_coef,a_opt_bias);

figure;
subplot(1,2,2)
plot(y,'Linewidth',1)
hold on
plot(y_hat_bias,'Linewidth',1)
grid on
grid minor
title('LMS prediction using scaled tanh activation function','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]','biased-tanh LMS estimate','Interpreter','latex')
ylim([-50 50])

subplot(1,2,1)
plot(y,'Linewidth',1)
hold on
plot(y_hat_bias,'Linewidth',1)
grid on
grid minor
title('First 200 time-steps','Fontsize',11)
xlabel('Sample index','Fontsize',11)
ylabel('Magnitude','Fontsize',11)
legend('y[n]','biased-tanh LMS estimate','Interpreter','latex')
ylim([-50 50])
xlim([0 200])

figure;
plot((error_bias).^2)
%Mean Square Error (dB)
MSE_opt_bias = pow2db(mean(abs(error_bias).^2)); %In dB
%Prediction Gain (dB)
Rp_opt_bias = pow2db(var(y_hat_bias)/var(error_bias));

figure;
subplot(1,2,1)
plot(a,MSE_bias)
subplot(1,2,2)
plot(a,Rp_bias)

%% Percentage difference between bias and no-bias 
Perc_diff_MSE=(abs(MSE_opt_bias - MSE_opt)/MSE_opt) *100
Perc_diff_Rp=(abs(Rp_opt_bias - Rp_opt)/Rp_opt) *100

figure; 
plot(y_hat_bias)
hold on
plot(y_hat_opt)

figure;
plot(weights_bias(2,:))
hold on
plot(weights_opt(1,:))

%% Finding MSE and Rp prior and post convergence
MSE_prior_bias = pow2db(mean(abs(error_bias(1:200)).^2)); %In dB
Rp_prior_bias = pow2db(var(y_hat_bias(1:200))/var(error_bias(1:200)));

MSE_post_bias = pow2db(mean(abs(error_bias(201:end)).^2)); %In dB
Rp_post_bias = pow2db(var(y_hat_bias(201:end))/var(error_bias(201:end)));


%% Observing convergence
figure;
for i=1:size(weights_bias,1)
    plot(weights_bias(i,:),'Linewidth',1)
    hold on
end
xline(150,'k--','Linewidth',1);
xlabel('Sample index','FontSize',11)
ylabel('Weight amplitude','FontSize',11)
title('Weight evolution for non-zero mean data','FontSize',11)
ylim([-0.025 0.025])
labels_2 = (0:size(weights_bias,1)-1);
legend(strcat('w_{',num2str(labels_2'),'}'),'Orientation','horizontal','Location','best')
grid on 
grid minor