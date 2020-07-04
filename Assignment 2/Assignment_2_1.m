%% Assignment 2.1
close all
clear all
clc

%% ==================================================================================
%% Part b)
%% ==================================================================================

% Synthesis Stage - synthesis of signal from white noise using AR model
N=1000;
noise_var = 0.25;
%Generating WGN with 0.25 variance
h=sqrt(noise_var).*randn(N,1); 

%Constructing AR model - Original sequence x(n)
% We want x(n)-a1*x(n-1)-a2*x(n-2)=h(n)
a=[1 -0.1 -0.8];
b=1;
x=filter(b,a,h); %Output of filter

% Analysis Stage - adaptation of coefficients to recreate original signal

%% One realisation

num_coef=length(a)-1;
mu1=0.05;
mu2=0.005;

% LMS predictor
[x_hat1, error1,weights1] = adaptive_lms(x,mu1,num_coef);
[x_hat2, error2,weights2] = adaptive_lms(x,mu2,num_coef);

figure(1);
subplot(1,2,1)
plot(pow2db(error1.^2),'Linewidth',1)
hold on
plot(pow2db(error2.^2),'Linewidth',1)
grid on
xlabel('Sample','FontSize',11); ylabel('Squared error (dB)','FontSize',11)
legend('\mu=0.05','\mu=0.01','FontSize',9)
title('1 realisation of x(n) used','FontSize',11)

%% Different realisations of x(n)
Num_iter = 100;
num_coef=length(a)-1;
mu1=0.05;
mu2=0.01;
%Initialising vectors/matrices
error_05=zeros(Num_iter,N);
error_01=zeros(Num_iter,N);
coefficients_05_1=zeros(Num_iter,N);
coefficients_05_2=zeros(Num_iter,N);
coefficients_01_1=zeros(Num_iter,N);
coefficients_01_2=zeros(Num_iter,N);

for iter=1:Num_iter
    
    % Synthesis Stage - synthesis of signal from white noise using AR model
    %Generating WGN with 0.25 variance
    h=sqrt(0.25).*randn(N,1); 
    %Constructing AR model - Original sequence x(n)
    a=[1 -0.1 -0.8];
    b=1;
    x=filter(b,a,h); %Output of filter

    % Analysis Stage - adaptation of coefficients to recreate original signal
    
    % LMS predictor
    [x_hat1, error1,weights1] = adaptive_lms(x,mu1,num_coef);
    [x_hat2, error2,weights2] = adaptive_lms(x,mu2,num_coef);
    
    %Storing the errors in a matrix
    error_05(iter,:)=error1;
    error_01(iter,:)=error2;
       
    %Storing the coefficients in a matrix
    %Negative sign needed due to inversion by adaptive_lms algorithm
    coefficients_05_1(iter,:) = weights1(1,:);
    coefficients_05_2(iter,:) = weights1(2,:);
    coefficients_01_1(iter,:) = weights2(1,:);
    coefficients_01_2(iter,:) = weights2(2,:);
end

%Obtaining the error squared in dB
error_05_mean=pow2db(mean(abs(error_05).^2));
error_01_mean=pow2db(mean(abs(error_01).^2));

figure(1);
subplot(1,2,2)
plot(error_05_mean,'Linewidth',1)
hold on
plot(error_01_mean,'Linewidth',1)
grid on
xlabel('Sample','FontSize',11); ylabel('Squared error (dB)','FontSize',11)
legend('\mu=0.05','\mu=0.01','FontSize',9)
title('100 realisations of x(n) used','FontSize',11)


%% ==================================================================================
%% Part c)
%% ==================================================================================

t_steady_05 = 500;% Time where transient response ends
t_steady_01 = 700;% Time where transient response ends
%Steady-state Learning curves
MSE_05 = mean(error_05(:,(t_steady_05:end)).^2);
MSE_01 = mean(error_01(:,(t_steady_01:end)).^2);

EMSE_05 = mean(MSE_05)-noise_var;
% % EMSE_05  = db2pow(EMSE_05);
EMSE_01 = mean(MSE_01)-noise_var;
% % EMSE_01  = db2pow(EMSE_01);
% % 
M_05 = EMSE_05/noise_var;
M_01 = EMSE_01/noise_var;

% Obtaining theoretical value 
R=[25/27 25/54; 25/54 25/27];
M_05_true = (mu1/2)*trace(R);
M_01_true = (mu2/2)*trace(R);


%% ==================================================================================
%% Part d)
%% ==================================================================================

%Obtaining the ensemble average of coefficients
a_pred_05_1= mean(coefficients_05_1);
a_pred_05_2= mean(coefficients_05_2);

a_pred_01_1= mean(coefficients_01_1);
a_pred_01_2= mean(coefficients_01_2);

%Plotting time evolution of coefficients
figure;
subplot(1,2,1)
plot(a_pred_05_1,'Linewidth',1.2)
hold on
yline(0.1,'--','color',[0 0.6 0],'Linewidth',1.2);
plot(a_pred_05_2,'Linewidth',1.2)
yline(0.8,'k--','Linewidth',1.2);
xlabel('Sample')
ylabel('Magnitude')
ylim([0 1])
title('Time evolution of Filter coefficients - \mu=0.05')
legend('$\hat{a}_{1}$','$a_{1}$',...
       '$\hat{a}_{2}$','$a_{2}$','Interpreter','latex');

subplot(1,2,2)
plot(a_pred_01_1,'Linewidth',1.2)
hold on
yline(0.1,'--','color',[0 0.6 0],'Linewidth',1.2);
plot(a_pred_01_2,'Linewidth',1.2)
yline(0.8,'k--','Linewidth',1.2);
xlabel('Sample')
ylabel('Magnitude')
ylim([0 1])
title('Time evolution of Filter coefficients - \mu=0.01')
legend('$\hat{a}_{1}$','$a_{1}$',...
       '$\hat{a}_{2}$','$a_{2}$','Interpreter','latex');
   
% Steady-state values of coefficients
a_pred_05_1_mean = mean(a_pred_05_1(t_steady_05:end));
a_pred_05_2_mean = mean(a_pred_05_2(t_steady_05:end));
a_pred_01_1_mean = mean(a_pred_01_1(t_steady_01:end));
a_pred_01_2_mean = mean(a_pred_01_2(t_steady_01:end));

%Considering bias and variance
a_pred_05_1_var= var(coefficients_05_1);
a_pred_01_1_var= var(coefficients_01_1);
a_pred_05_2_var= var(coefficients_05_2);
a_pred_01_2_var= var(coefficients_01_2);

a_pred_05_1_bias = mean(coefficients_05_1) - (-a(2));
a_pred_01_1_bias = mean(coefficients_01_1) - (-a(2));
a_pred_05_2_bias = mean(coefficients_05_2) - (-a(3));
a_pred_01_2_bias = mean(coefficients_01_2) - (-a(3));

figure;
subplot(2,2,1)
plot(a_pred_05_1_var,'Linewidth',1)
hold on
plot(a_pred_01_1_var,'Linewidth',1)
grid on
xlabel('Sample')
ylabel('Variance')
legend('\mu=0.05','\mu=0.01')
title('Variance of $\hat{a}_{1}$','Interpreter','latex','FontSize',11)

subplot(2,2,2)
plot(a_pred_05_2_var,'Linewidth',1)
hold on
plot(a_pred_01_2_var,'Linewidth',1)
grid on
xlabel('Sample')
ylabel('Variance')
legend('\mu=0.05','\mu=0.01')
title('Variance of $\hat{a}_{2}$','Interpreter','latex','FontSize',11)

subplot(2,2,3)
plot(abs(a_pred_05_1_bias),'Linewidth',1)
hold on
plot(abs(a_pred_01_1_bias),'Linewidth',1)
grid on
xlabel('Sample')
ylabel('Bias')
legend('\mu=0.05','\mu=0.01')
title('Bias of $\hat{a}_{1}$','Interpreter','latex','FontSize',11)

subplot(2,2,4)
plot(abs(a_pred_05_2_bias),'Linewidth',1)
hold on
plot(abs(a_pred_01_2_bias),'Linewidth',1)
grid on
xlabel('Sample')
ylabel('Bias')
legend('\mu=0.05','\mu=0.01')
title('Bias of $\hat{a}_{1}$','Interpreter','latex','FontSize',11)

%Plotting MSE
a_pred_05_1_mse= (a_pred_05_1_bias.^2)+a_pred_05_1_var;
a_pred_01_1_mse= (a_pred_01_1_bias.^2)+a_pred_01_1_var;
a_pred_05_2_mse= (a_pred_05_2_bias.^2)+a_pred_05_2_var;
a_pred_01_2_mse= (a_pred_01_2_bias.^2)+a_pred_01_2_var;

figure;
subplot(1,2,1)
plot(a_pred_05_1_mse,'Linewidth',1)
hold on
plot(a_pred_01_1_mse,'Linewidth',1)
grid on
xlabel('Sample')
ylabel('MSE')
legend('\mu=0.05','\mu=0.01')
title('MSE of LMS - a_{1}')

subplot(1,2,2)
plot(a_pred_05_2_mse,'Linewidth',1)
hold on
plot(a_pred_01_2_mse,'Linewidth',1)
grid on
xlabel('Sample')
ylabel('MSE')
legend('\mu=0.05','\mu=0.01')
title('MSE of LMS - a_{2}')


%% ==================================================================================
%% Part f)
%% ==================================================================================

%% Different realisations of x(n)
Num_iter = 100;

%AR true model coefficients
a=[1 -0.1 -0.8];
b=1;

mu=[0.01 0.05];
gamma = [0.2 0.5 0.8];
%Matrix for estimated coefficients
coefficients_leaky=zeros(length(a)-1,Num_iter,N);

for iter=1:Num_iter
    
    %selecting mu and gamma
    m=mu(2);
    g=0; %gamma(3);
    % Synthesis Stage - synthesis of signal from white noise using AR model
    h=sqrt(0.25).*randn(N,1); %This is already standardised
    x=filter(b,a,h); %Output of filter

    % Analysis Stage - adaptation of coefficients to recreate original signal

    num_coef=length(a)-1;
        
    % leaky-LMS predictor
    [~,~,weights] = leaky_lms(x,m,g,num_coef);
         
    %Storing the coefficients in a matrix
    %Negative sign needed due to inversion by adaptive_lms algorithm
    coefficients_leaky(1,iter,:) = weights(1,:);
    coefficients_leaky(2,iter,:) = weights(2,:);
end

%Obtaining the ensemble average of coefficients
a1_leaky(:)= mean(coefficients_leaky(1,:,:));
a2_leaky(:)= mean(coefficients_leaky(2,:,:));

% Steady-state values of coefficients
t_steady = 600;
a1_leaky_mean = mean(a1_leaky(t_steady:end))
a2_leaky_mean = mean(a2_leaky(t_steady:end))

absolute_error_1=(abs(a1_leaky_mean-(-a(2)))/(-a(2)))*100
absolute_error_2=(abs(a2_leaky_mean-(-a(3)))/(-a(3)))*100

figure;
plot(a1_leaky,'Linewidth',1.2)
hold on
yline(0.1,'--','color',[0 0.6 0],'Linewidth',1.2);
plot(a2_leaky,'Linewidth',1.2)
yline(0.8,'k--','Linewidth',1.2);
xlabel('Sample')
ylabel('Magnitude')
ylim([0 1])
title(['Time evolution of Filter coefficients - \mu=',num2str(m),' and \gamma=',num2str(g)])
legend('$\hat{a}_{1}$','$a_{1}$',...
       '$\hat{a}_{2}$','$a_{2}$','Interpreter','latex');
   
figure;
x_hat = [];
x_hat(1)=0;
x_hat(2)=0;
for n=3:N
    x_hat(n) = a1_leaky_mean*x(n-1) + a2_leaky_mean*x(n-2);
end
plot(x)
hold on
plot(x_hat)

%% Estimating optimal solutions
Num_iter = 1000;
N=1000;

%AR true model coefficients
a=[1 -0.1 -0.8];
b=1;

for iter=1:Num_iter
    % Synthesis Stage - synthesis of signal from white noise using AR model
    h=sqrt(0.25).*randn(N,1); %This is already standardised
    x=filter(b,a,h); %Output of filter
    
    for n = 2 : N-1
        x_n = x(n:-1:n-2+1); %Obtaining past values of input
        %Approximation to cross-correlation vector
        p(:,n) = x(n+1)*x_n;
    end
    p_steady = mean(p(:,600:end)')';
    p_steady_mtx(:,iter) = p_steady;
end

p_steady_mean = mean(p_steady_mtx')';
I=eye(2,2);
gamma=0.99;
R_inverse = inv(R+gamma.*I);
wopt = R_inverse*p_steady_mean

%% Considering weight misalignment
a1_leaky_bias(:) = mean(coefficients_leaky(1,:,:))-0.1;
a1_leaky_variance(:)= var(coefficients_leaky(1,:,:));

a2_leaky_bias(:) = mean(coefficients_leaky(2,:,:))-0.8;
a2_leaky_variance(:)= var(coefficients_leaky(2,:,:));

a1_leaky_mse= (a1_leaky_bias.^2)+a1_leaky_variance;
a2_leaky_mse= (a2_leaky_bias.^2)+a2_leaky_variance;

figure;
subplot(1,2,1)
plot(a1_leaky_mse,'Linewidth',1)
xlabel('Sample')
ylabel('MSE')
title('MSE for coefficient 1');

subplot(1,2,2)
plot(a2_leaky_mse,'Linewidth',1)
xlabel('Sample')
ylabel('MSE')
title('MSE for coefficient 2');

mean(a1_leaky_mse(600:end))
mean(a2_leaky_mse(600:end))