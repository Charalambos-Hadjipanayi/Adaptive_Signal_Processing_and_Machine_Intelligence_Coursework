%% Assignment 3.1 e
close all; clear all ; clc

% Time axis (n)
N=5000; %time samples to be used
n=0:N-1;

% System frequency
f_0 = 50; %Hz
% Sampling frequency
f_s = 1000; %Hz

%% Generating balanced system - All voltages the same and delta=0.
V_bal=[1,1,1]; %Va,Vb,Vc
delta_bal =[0,0]; %delta_b, delta_c
phi_bal =0; %Nominal phase shift-common for all voltages

% The phase voltages
V_a_bal = V_bal(1)*cos(2*pi*(f_0/f_s)*n + phi_bal);
V_b_bal = V_bal(2)*cos(2*pi*(f_0/f_s)*n + phi_bal+delta_bal(1) - (2*pi/3));
V_c_bal = V_bal(3)*cos(2*pi*(f_0/f_s)*n + phi_bal+delta_bal(2) + (2*pi/3));

%Applying Clarke Transform
Voltage_TF_bal = clarke_tf([V_a_bal ; V_b_bal ; V_c_bal]);

%Forming Clarke Voltage
Voltage_Clarke_bal = Voltage_TF_bal(2,:) + 1j*Voltage_TF_bal(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_bal= mean(Voltage_Clarke_bal .^ 2);
covariance_bal = mean(abs(Voltage_Clarke_bal) .^ 2);
circularity_coef_bal = abs(pseudocovariance_bal)/covariance_bal;

%% Generating unbalanced system - Voltages vary and delta varies.
V_unbal=[1,0.8,1.2]; %Va,Vb,Vc
delta_unbal =[0.4,0.9]; %delta_b, delta_c
phi_unbal =0; %Nominal phase shift-common for all voltages

% The phase voltages
V_a_unbal = V_unbal(1)*cos(2*pi*(f_0/f_s)*n + phi_unbal);
V_b_unbal = V_unbal(2)*cos(2*pi*(f_0/f_s)*n + phi_unbal+delta_unbal(1) - (2*pi/3));
V_c_unbal = V_unbal(3)*cos(2*pi*(f_0/f_s)*n + phi_unbal+delta_unbal(2) + (2*pi/3));

%Applying Clarke Transform
Voltage_TF_unbal = clarke_tf([V_a_unbal ; V_b_unbal ; V_c_unbal]);

%Forming Clarke Voltage
Voltage_Clarke_unbal = Voltage_TF_unbal(2,:) + 1j*Voltage_TF_unbal(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_unbal= mean(Voltage_Clarke_unbal .^ 2);
covariance_unbal = mean(abs(Voltage_Clarke_unbal) .^ 2);
circularity_coef_unbal = abs(pseudocovariance_unbal)/covariance_unbal;

%% Part e) - CLMS and ACLMS algorithms to detect frequency
coef=1; %AR(1) process
mu=0.05;

%Balanced system
[~,error_bal_CLMS,h_bal_CLMS] = CLMS_pred_AR(Voltage_Clarke_bal,mu,coef); 
f0_bal_CLMS = abs((f_s/(2*pi))*atan(imag(h_bal_CLMS)./real(h_bal_CLMS)));

% [~,error_bal_ACLMS,h_bal_ACLMS,g_bal_ACLMS] = aclms_ar(Voltage_Clarke_bal,mu,coef); 
% f0_bal_ACLMS = f_s/(2*pi)*atan( sqrt( imag(h_bal_ACLMS).^2 - abs(g_bal_ACLMS).^2 ) ./ real(h_bal_ACLMS) );

[~,error_bal_ACLMS,h_balanced_aclms,g_balanced_aclms]=ACLMS_pred_AR(Voltage_Clarke_bal,mu,coef);
f0_bal_ACLMS=abs((f_s/(2*pi))*atan((sqrt( (imag(h_balanced_aclms).^2) - (abs(g_balanced_aclms).^2) ))./ (real(h_balanced_aclms))));

figure;
subplot(1,2,1)
plot(f0_bal_CLMS,'Linewidth',1.2)
hold on
plot(f0_bal_ACLMS,'Linewidth',1.2)
yline(50,'k--','Linewidth',1);
xlabel('Sample index (n)','FontSize',11)
ylabel('Estimated frequency f_{0}','FontSize',11)
title('Balanced system - Frequency estimation','FontSize',11)
legend('CLMS','ACLMS','Theoretical f_{0}')
ylim([0 60])
xlim([0 500])
grid on
grid minor

%% Unbalanced system - both phase and magnitude
[~,error_unbal_CLMS,h_unbal_CLMS] = CLMS_pred_AR(Voltage_Clarke_unbal,mu,coef); 
f0_unbal_CLMS = abs(f_s/(2*pi)*atan(imag(h_unbal_CLMS)./real(h_unbal_CLMS)));

[~,error_unbal_ACLMS,h_unbal_ACLMS,g_unbal_ACLMS] = ACLMS_pred_AR(Voltage_Clarke_unbal,mu,coef); 
f0_unbal_ACLMS = abs(f_s/(2*pi)*atan(sqrt((imag(h_unbal_ACLMS)).^2 - abs(g_unbal_ACLMS).^2)./real(h_unbal_ACLMS)));

subplot(1,2,2)
plot(f0_unbal_CLMS,'Linewidth',1)
hold on
plot(f0_unbal_ACLMS,'Linewidth',1)
yline(50,'k--','Linewidth',1);
xlabel('Sample index (n)','FontSize',11)
ylabel('Estimated frequency f_{0}','FontSize',11)
title('Unbalanced system - Frequency estimation','FontSize',11)
legend('CLMS','ACLMS','Theoretical f_{0}')
ylim([0 60])
xlim([0 500])
grid on
grid minor

%% Output frequencies
f0_bal_clms = mean(f0_bal_CLMS(N-200:end))
f0_bal_aclms = mean(f0_bal_ACLMS(N-200:end))
f0_unbal_clms = mean(f0_unbal_CLMS(N-200:end))
f0_unbal_aclms = mean(f0_unbal_ACLMS(N-200:end))


