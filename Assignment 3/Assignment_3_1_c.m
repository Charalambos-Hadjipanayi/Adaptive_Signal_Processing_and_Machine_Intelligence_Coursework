%% Assignment 3.1 c
close all; clear all ; clc

%% Generating balanced system - All voltages the same and delta=0.
V=[2,2,2]; %Va,Vb,Vc
delta =[0,0]; %delta_b, delta_c
phi =0; %Nominal phase shift-common for all voltages

% Time axis (n)
N=100000; %time samples to be used
n=0:N-1;

% System frequency
f_0 = 50; %Hz
% Sampling frequency
f_s = 10000; %Hz

% The phase voltages
V_a = V(1)*cos(2*pi*(f_0/f_s)*n + phi);
V_b = V(2)*cos(2*pi*(f_0/f_s)*n + phi+delta(1) - (2*pi/3));
V_c = V(3)*cos(2*pi*(f_0/f_s)*n + phi+delta(2) + (2*pi/3));

%Forming vector of voltages
V_triphase = [V_a ; V_b ; V_c];

%Applying Clarke Transform
Voltage_TF = clarke_tf(V_triphase);

%Forming Clarke Voltage
Voltage_Clarke = Voltage_TF(2,:) + 1j*Voltage_TF(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_bal= mean(Voltage_Clarke .^ 2);
covariance_bal = mean(abs(Voltage_Clarke) .^ 2);
circularity_coef_bal = abs(pseudocovariance_bal)/covariance_bal;

%Circularity plot
figure(1);
subplot(1,2,1)
scatter(real(Voltage_Clarke),imag(Voltage_Clarke),'.')
xlabel('Real part - $v_{\alpha}$','Interpreter','Latex','FontSize',10)
ylabel('Imaginary part - $v_{\beta}$','Interpreter','Latex','FontSize',10)
grid on
grid minor
title({
    ['Balanced system with |\rho| =',num2str(round(circularity_coef_bal,4))] 
    ['V_{a}=V_{b}=V_{c}=',num2str(V(3)),',\Delta_{b}=\Delta_{c}=',num2str(delta(2))] 
    });


%% Generating unbalanced system - Voltages vary and delta varies.
V=[1,3,5]; %Va,Vb,Vc
delta =[0.2,0.8]; %delta_b, delta_c
phi =0; %Nominal phase shift-common for all voltages

% Time axis (n)
N=100000; %time samples to be used
n=0:N-1;

% System frequency
f_0 = 50; %Hz
% Sampling frequency
f_s = 10000; %Hz

% The phase voltages
V_a = V(1)*cos(2*pi*(f_0/f_s)*n + phi);
V_b = V(2)*cos(2*pi*(f_0/f_s)*n + phi+delta(1) - (2*pi/3));
V_c = V(3)*cos(2*pi*(f_0/f_s)*n + phi+delta(2) + (2*pi/3));

%Forming vector of voltages
V_triphase = [V_a ; V_b ; V_c];

%Applying Clarke Transform
Voltage_TF = clarke_tf(V_triphase);

%Forming Clarke Voltage
Voltage_Clarke = Voltage_TF(2,:) + 1j*Voltage_TF(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_unbal= mean(Voltage_Clarke .^ 2);
covariance_unbal = mean(abs(Voltage_Clarke) .^ 2);
circularity_coef_unbal = abs(pseudocovariance_unbal)/covariance_unbal;

%Circularity plot
figure (1);
subplot(1,2,2)
scatter(real(Voltage_Clarke),imag(Voltage_Clarke),'r','.')
xlabel('Real part - $v_{\alpha}$','Interpreter','Latex','FontSize',10)
ylabel('Imaginary part - $v_{\beta}$','Interpreter','Latex','FontSize',10)
grid on
grid minor
title({
    ['Unbalanced Magnitude and phase with |\rho| =',num2str(round(circularity_coef_unbal,4))] 
    ['V_{a}=',num2str(V(1)),',V_{b}=',num2str(V(2)),',V_{c}=',num2str(V(3)),',\Delta_{b}=',num2str(delta(1)),',\Delta_{c}=',num2str(delta(2))] 
    });

%% Generating unbalanced system - Voltages vary and delta=0.
V=[1,0.5,0.8,1,0.5,1.4,1,2,0.6]; %Va,Vb,Vc
delta =[0,0]; %delta_b, delta_c
phi =0; %Nominal phase shift-common for all voltages

% Time axis (n)
N=100000; %time samples to be used
n=0:N-1;

% System frequency
f_0 = 50; %Hz
% Sampling frequency
f_s = 10000; %Hz

circularity_coef_unbal=zeros(length(V)/3,1);
N_iter = length(circularity_coef_unbal);

for i =0:N_iter-1
% The phase voltages
V_a = V(1+i)*cos(2*pi*(f_0/f_s)*n + phi);
V_b = V(N_iter+1+i)*cos(2*pi*(f_0/f_s)*n + phi+delta(1) - (2*pi/3));
V_c = V(2*(N_iter)+1+i)*cos(2*pi*(f_0/f_s)*n + phi+delta(2) + (2*pi/3));

%Forming vector of voltages
V_triphase = [V_a ; V_b ; V_c];

%Applying Clarke Transform
Voltage_TF = clarke_tf(V_triphase);

%Forming Clarke Voltage
Voltage_Clarke = Voltage_TF(2,:) + 1j*Voltage_TF(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_unbal= mean(Voltage_Clarke .^ 2);
covariance_unbal = mean(abs(Voltage_Clarke) .^ 2);
circularity_coef_unbal(i+1) = abs(pseudocovariance_unbal)/covariance_unbal;

%Circularity plot
figure (2);
subplot(1,2,1)
scatter(real(Voltage_Clarke),imag(Voltage_Clarke),'.')
hold on
end
xlabel('Real part - $v_{\alpha}$','Interpreter','Latex','FontSize',10)
ylabel('Imaginary part - $v_{\beta}$','Interpreter','Latex','FontSize',10)
grid on
grid minor
title(['Unbalanced Magnitude - \Delta_{b}=',num2str(delta(1)),',\Delta_{c}=',num2str(delta(2))])

%Forming Voltage magnitude used in column vectors
V_1 = V(1:N_iter)';
V_2 = V(N_iter+1:2*N_iter)';
V_3 = V(2*(N_iter)+1:3*N_iter)';

legend(strcat('V_{a}=',num2str(V_1),', V_{b}=',num2str(V_2),', V_{c}=',num2str(V_3),',|\rho|=',num2str(round(circularity_coef_unbal,3))),'Fontsize',8)
ylim([-3 3])
xlim([-3 3])

%% Generating unbalanced system - Same voltages vary but vary delta.
V=[1,1,1]; %Va,Vb,Vc
delta =[0,1.5,0.7,0,0.4,1.2]; %delta_b, delta_c
phi =0; %Nominal phase shift-common for all voltages

% Time axis (n)
N=100000; %time samples to be used
n=0:N-1;

% System frequency
f_0 = 50; %Hz
% Sampling frequency
f_s = 10000; %Hz

circularity_coef_unbal=zeros(length(delta)/2,1);
N_iter = length(circularity_coef_unbal);

for i =0:N_iter-1
% The phase voltages
V_a = V(1)*cos(2*pi*(f_0/f_s)*n + phi);
V_b = V(2)*cos(2*pi*(f_0/f_s)*n + phi+delta(1+i) - (2*pi/3));
V_c = V(3)*cos(2*pi*(f_0/f_s)*n + phi+delta(N_iter+1+i) + (2*pi/3));

%Forming vector of voltages
V_triphase = [V_a ; V_b ; V_c];

%Applying Clarke Transform
Voltage_TF = clarke_tf(V_triphase);

%Forming Clarke Voltage
Voltage_Clarke = Voltage_TF(2,:) + 1j*Voltage_TF(3,:);

% Finding the Circularity coefficient for balanced system
pseudocovariance_unbal= mean(Voltage_Clarke .^ 2);
covariance_unbal = mean(abs(Voltage_Clarke) .^ 2);
circularity_coef_unbal(i+1) = abs(pseudocovariance_unbal)/covariance_unbal;

%Circularity plot
figure (2);
subplot(1,2,2)
scatter(real(Voltage_Clarke),imag(Voltage_Clarke),'.')
hold on
end
xlabel('Real part - $v_{\alpha}$','Interpreter','Latex','FontSize',10)
ylabel('Imaginary part - $v_{\beta}$','Interpreter','Latex','FontSize',10)
grid on
grid minor
title(['Unbalanced Phase - V_{a}=V_{b}=V_{c}=',num2str(V(3))])

%Forming Voltage magnitude used in column vectors
delta_1 = delta(1:N_iter)';
delta_2 = delta(N_iter+1:end)';

legend(strcat('\Delta_{b}=',num2str(delta_1),', \Delta_{c}=',num2str(delta_2),',|\rho|=',num2str(round(circularity_coef_unbal,3))),'Fontsize',8)
ylim([-3 3])
xlim([-3 3])
