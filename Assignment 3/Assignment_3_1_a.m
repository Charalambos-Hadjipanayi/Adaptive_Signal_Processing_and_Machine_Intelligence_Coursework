%% Assignment 3.1 a
close all; clear all ; clc

%Signal length
N=100000;
%Generating circular white noise
x = randn(1,N)+1j*randn(1,N); 

%Generating WLMA(1) process
signal = zeros(1,length(x));
b = [1 (1.5+1*1j) (2.5-0.5*1j)];
signal(1) = x(1);
for i = 1:N-1
    signal(i+1) = b(1)*x(i+1) + b(2)*x(i) + b(3)*conj(x(i));
end

% Finding the Circularity coefficient for white noise
pseudocovariance_WGN= mean(x .^ 2);
covariance_WGN = mean(abs(x) .^ 2);
circularity_coef_WGN = abs(pseudocovariance_WGN)/covariance_WGN;

% Finding the Circularity coefficient for WLMA(1) process
pseudocovariance_signal= mean(signal .^ 2);
covariance_signal = mean(abs(signal) .^ 2);
circularity_coef_signal = abs(pseudocovariance_signal)/covariance_signal;

%Checking circularity of WGN
figure;
subplot(1,2,1)
scatter(real(x),imag(x),'.')
xlabel('Real part')
ylabel('Imaginary part')
grid on
grid minor
title(['Circular WGN with |\rho| =',num2str(round(circularity_coef_WGN,4))]) 

%Checking circularity of Signal (WLMA(1) process)
subplot(1,2,2)
scatter(real(signal),imag(signal),'k','.')
xlabel('Real part')
ylabel('Imaginary part')
grid on
grid minor
title(['WLMA(1) with |\rho| =',num2str(round(circularity_coef_signal,3))]) 


%% Considering 100 iterations
N_iter=100;
mu = 0.03;  %learning rate
num_coef=length(b)-1;
N=1000; %Signal length

%Initialising matrices
error_CLMS = zeros(N_iter,N);
error_ACLMS = zeros(N_iter,N);

for iter=1:N_iter
    %Generating circular white noise
    x = randn(1,N)+1j*randn(1,N);
    %Generating WLMA(1) process
    signal = zeros(1,length(x));
    b = [1 (1.5+1*1j) (2.5-0.5*1j)];
    signal(1) = x(1);
    for i = 1:N-1
        signal(i+1) = b(1)*x(i+1) + b(2)*x(i) + b(3)*conj(x(i));
    end
        
    % Implementing the CLMS algorithm
    [signal_est_CLMS,error_1,~] = CLMS_MA(x,signal,mu,num_coef);
    error_CLMS(iter,:) = error_1;
    
    %Implementing the ACLMS algorithm
    [signal_est_ACLMS,error_2,~,~] = ACLMS_MA(x,signal,mu,num_coef);
    error_ACLMS(iter,:) = error_2;
end

%Obtaining learning curves
figure;
subplot(1,3,1)
plot(pow2db(mean(abs(error_CLMS).^2,1)),'Linewidth',1)
hold on
plot(pow2db(mean(abs(error_ACLMS).^2,1)),'Linewidth',1)
xlabel('Sample Index','FontSize',11);
ylabel('Square error (dB)','FontSize',11);
legend('CLMS','ACLMS')
title('Learning curves for CLMS and ACLMS algorithms','FontSize',11)
grid on
grid minor

subplot(1,3,2)
scatter(real(signal),imag(signal),'.')
hold on
scatter(real(signal_est_CLMS),imag(signal_est_CLMS),'.')
xlabel('Real part','FontSize',11)
ylabel('Imaginary part','FontSize',11)
title('One realisation - CLMS algorithm','FontSize',11)
legend('Original process','CLMS estimate', 'Orientation','horizontal')
grid on
grid minor

subplot(1,3,3)
scatter(real(signal),imag(signal),'.')
hold on
scatter(real(signal_est_ACLMS),imag(signal_est_ACLMS),'.')
xlabel('Real part','FontSize',11)
ylabel('Imaginary part','FontSize',11)
title('One realisation - ACLMS algorithm','FontSize',11)
legend('Original process','ACLMS estimate', 'Orientation','horizontal')
grid on
grid minor

    