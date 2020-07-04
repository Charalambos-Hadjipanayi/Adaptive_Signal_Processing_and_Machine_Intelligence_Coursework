%% Assignment 2.2
close all; clear all ; clc

%% ==================================================================================
%% Part a)
%% ==================================================================================

Num_iter = 100; %number of iterations
N=2000; %Data samples

%Model Parameters
a=[1];
b=[1 0.9];
num_coef=length(b); %Filter length
noise_var=0.5;
mu_0=0.1; %Initial step-size for GASS
mu=[0.01,0.1];
rho=0.005;
alpha = 0.9;

%Initialising matrices to store error
error_1=zeros(N,Num_iter);
error_2=zeros(N,Num_iter);
error_3=zeros(N,Num_iter);
error_4=zeros(N,Num_iter);
error_5=zeros(N,Num_iter);

%Initialising matrices to store weight wo
wmtx_1=zeros(num_coef-1,N);
wmtx_2=zeros(num_coef-1,N);
wmtx_3=zeros(num_coef-1,N);
wmtx_4=zeros(num_coef-1,N);
wmtx_5=zeros(num_coef-1,N);

for iter=1:100
    % Generating one signal only
    h = random('Normal', 0, noise_var, N, 1);
    %Constructing MA model - Original sequence x(n)
    % We want x(n)=a1*h(n-1)+h(n)
    x=filter(b,a,h); %Output of filter
    
    % 1. Standard LMS with mu=0.01
    [~,er_1,weights_1] = adaptive_lms_MA(h,x,mu(1),num_coef);
    error_1(:,iter)=er_1;
    wmtx_1(iter,:)=weights_1(1,:);
        
    % 2. Standard LMS with mu=0.1
    [~,er_2,weights_2] = adaptive_lms_MA(h,x,mu(2),num_coef);
    error_2(:,iter)=er_2;
    wmtx_2(iter,:)=weights_2(1,:);
    
    % 3. Benveniste Algorithm
    [~,~,er_3,weights_3] = benveniste(h,x,mu_0, rho,num_coef);
    error_3(:,iter)=er_3;
    wmtx_3(iter,:)=weights_3(1,:);
    
    % 4. Ang & Farhang Algorithm
    [~,~,er_4,weights_4] = ang_farhang(h,x,mu_0, rho,alpha,num_coef);
    error_4(:,iter)=er_4;
    wmtx_4(iter,:)=weights_4(1,:);
    
    % 5. Matthews % Xie Algorithm
    [~,~,er_5,weights_5] = matthews_xie(h,x,mu_0, rho,num_coef);
    error_5(:,iter)=er_5;
    wmtx_5(iter,:)=weights_5(1,:);    
end

%Obtaining weights errors
w0=0.9*ones(Num_iter,N); %Theoretical value
w_error_mtx_0_1 = w0 - wmtx_1;
w_error_mtx_0_2 = w0 - wmtx_2;
w_error_mtx_0_3 = w0 - wmtx_3;
w_error_mtx_0_4 = w0 - wmtx_4;
w_error_mtx_0_5 = w0 - wmtx_5;

figure;
plot(mean(wmtx_1,1),'Linewidth',1)
hold on
plot(mean(wmtx_2,1),'Linewidth',1)
plot(mean(wmtx_3,1),'Linewidth',1)
plot(mean(wmtx_4,1),'Linewidth',1)
plot(mean(wmtx_5,1),'Linewidth',1)
yline(0.9,'k--','Linewidth',1);
xlabel('Sample','FontSize',11)
ylabel('Magnitude','FontSize',11)
legend('\mu=0.01','\mu=0.1','Benveniste','Ang&Farhang','Matthews&Xie','Theoretical')
title('Time evolution of weight w_{0}','FontSize',11)

%% Obtaining weight error curves
figure(2);
subplot(1,3,1)
plot(mean(w_error_mtx_0_1,1),'Linewidth',0.8)
hold on
plot(mean(w_error_mtx_0_2,1),'Linewidth',0.8)
plot(mean(w_error_mtx_0_3,1),'Linewidth',0.8)
plot(mean(w_error_mtx_0_4,1),'Linewidth',0.8)
plot(mean(w_error_mtx_0_5,1),'Linewidth',0.8)
xlabel('Sample','FontSize',11)
ylabel('Weight error','FontSize',11)
legend('\mu=0.01','\mu=0.1','Benveniste','Ang&Farhang','Matthews&Xie')
title('Weight error curves of algorithms','FontSize',11)
grid on; grid minor
ylim([-0.1 0.9])

%% Obtaining weight error curves in dB
figure(2);
subplot(1,3,2)
plot(pow2db(mean(w_error_mtx_0_1.^2,1)),'Linewidth',1)
hold on
plot(pow2db(mean(w_error_mtx_0_2.^2,1)),'Linewidth',1)
plot(pow2db(mean(w_error_mtx_0_3.^2,1)),'Linewidth',1)
plot(pow2db(mean(w_error_mtx_0_4.^2,1)),'Linewidth',1)
plot(pow2db(mean(w_error_mtx_0_5.^2,1)),'Linewidth',1)
xlabel('Sample','FontSize',11)
ylabel('Squared weight error (dB)','FontSize',11)
legend('\mu=0.01','\mu=0.1','Benveniste','Ang&Farhang','Matthews&Xie')
title('Weight error curves of algorithms','FontSize',11)
grid on; grid minor

%% Loading method 1 data
load('method_1.mat');

figure(2);
subplot(1,3,3)
plot(method_1(1,:),'Linewidth',1)
hold on
plot(method_1(2,:),'Linewidth',1)
plot(method_1(3,:),'Linewidth',1)
xlabel('Sample','FontSize',11)
ylabel('Squared weight error (dB)','FontSize',11)
legend('\mu=0.01','\mu=0.1','\mu=1')
title('Method 1 - Idealised case','FontSize',11)
grid on; grid minor

%% Obtaining square error curves
square_error_1 = pow2db(mean(error_1.^2,2));
square_error_2 = pow2db(mean(error_2.^2,2));
square_error_3 = pow2db(mean(error_3.^2,2));
square_error_4 = pow2db(mean(error_4.^2,2));
square_error_5 = pow2db(mean(error_5.^2,2));

figure(3);
plot(square_error_1,'Linewidth',1.5)
hold on
plot(square_error_2,'Linewidth',1.5)
plot(square_error_3,'Linewidth',1.5)
plot(square_error_4,'Linewidth',1.5)
plot(square_error_5,'Linewidth',1.5)
xlabel('Sample','FontSize',11)
ylabel('Square Error (dB)','FontSize',11)
legend('\mu=0.01','\mu=0.1','Benveniste','Ang&Farhang','Matthews&Xie')
title('Square error of algorithms','FontSize',11)
grid on; grid minor

%% ==================================================================================
%% Part c)
%% ==================================================================================
Num_iter = 100; %number of iterations
N=1000; %Data samples

%Model Parameters
a=[1];
b=[1 0.9];
num_coef=length(b); %Filter length
noise_var=0.5;
mu_0=0.1; %Initial step-size for GASS
mu=[0.1 0.5]; %Initial step-size for GNGD
rho=0.01;

%Initialising matrices to store error
error_benveniste=zeros(N,Num_iter);
error_GNGD=zeros(N,Num_iter);

%Initialising matrices to store weight wo
wmtx_benveniste=zeros(num_coef-1,N);
wmtx_GNGD=zeros(num_coef-1,N);

for i=1:length(mu)
    epsilon_0=1/mu(i);
    for iter=1:100
        % Generating one signal only
        h = random('Normal', 0, noise_var, N, 1);
        %Constructing MA model - Original sequence x(n)
        % We want x(n)=a1*h(n-1)+h(n)
        x=filter(b,a,h); %Output of filter
        
        % 1. Benveniste algorithm
        [~,~,er_benveniste,weights_benveniste] = benveniste(h,x,mu_0,rho,num_coef);
        error_benveniste(:,iter)=er_benveniste;
        wmtx_benveniste(iter,:)=weights_benveniste(1,:);
        
        % 2. GNGD algorithm
        [~,er_GNGD,weights_GNGD] = GNGD(h,x,mu(i),epsilon_0,rho,num_coef);
        error_GNGD(:,iter)=er_GNGD;
        wmtx_GNGD(iter,:)=weights_GNGD(1,:);
    end
    
    %Obtaining mean weights - Ensemble-averaged weight evolutions
    w0=0.9*ones(1,N); %Theoretical value
    w0_benveniste = mean(wmtx_benveniste,1);
    w0_GNGD = mean(wmtx_GNGD,1);
    
    %Obtaining evolution of weight estimates
    figure(4)
    subplot(1,2,i)
    plot(w0_GNGD,'Linewidth',1.2)
    hold on
    plot(w0_benveniste,'Linewidth',1.2)
    yline(0.9,'k--','Linewidth',0.6);
    xlabel('Sample','FontSize',11)
    ylabel('Magnitude','FontSize',11)
    legend('GNGD','Benveniste','Theoretical')
    title(['Time evolution of weight w_{0} - \mu=' num2str(mu(i))],'FontSize',11)
    grid on; grid minor
    ylim([0 0.95]);
    
    %Obtaining weight error curves
    figure(5)
    subplot(1,2,i)
    plot(w0-w0_benveniste,'Linewidth',1.8)
    hold on
    plot(w0-w0_GNGD,'Linewidth',1.8)
    xlabel('Sample','FontSize',11)
    ylabel('Weight error','FontSize',11)
    legend('Benveniste','GNGD')
    title('Weight error curves of algorithms','FontSize',11)
    grid on; grid minor
    ylim([0 0.9]);
end

%% Checking performance with changing parameters - rho
N_iter=100;
rho=(0.0001:0.0001:0.01);
w0=0.9*ones(N_iter,N); %Theoretical value
mu_0=0.1; %Initial step-size for GASS
mu=0.1; %Initial step-size for GNGD
for rho_index =1:length(rho)
    rho_index
    for iter=1:N_iter
        
        % Generating one signal only
        h = random('Normal', 0, noise_var, N, 1);
        %Constructing MA model - Original sequence x(n)
        % We want x(n)=a1*h(n-1)+h(n)
        x=filter(b,a,h); %Output of filter
        
        % 1. Benveniste algorithm
        [~,~,er_benveniste,weights_benveniste] = benveniste_Method_1(h,x,mu_0,rho(rho_index),num_coef);
        error_benveniste(:,iter)=er_benveniste;
        wmtx_benveniste(iter,:)=weights_benveniste(1,:);
        
        % 2. GNGD algorithm
        [~,er_GNGD,weights_GNGD] = GNGD_Method_1(h,x,mu,epsilon_0,rho(rho_index),num_coef);
        error_GNGD(:,iter)=er_GNGD;
        wmtx_GNGD(iter,:)=weights_GNGD(1,:);
    end
    
    w_error_benveniste = mean(abs(w0 - wmtx_benveniste),1);
    w_error_GNGD = mean(abs(w0 - wmtx_GNGD),1);
    
    % Finding concergence time
    threshold = 0.001; %Setting threshold for error
    c_s_ben = find(w_error_benveniste<threshold);
    if isempty(c_s_ben)
        c_s_ben=N;
    end
    convergence_time_ben (rho_index)= min (c_s_ben);
    c_s_GNGD = find(w_error_GNGD<threshold);
    if isempty(c_s_GNGD)
        c_s_GNGD=N;
    end
    convergence_time_GNGD (rho_index)= min (c_s_GNGD);
end

%%
figure(6);
subplot(1,3,1)
plot(rho,1./convergence_time_ben,'Linewidth',1)
hold on
plot(rho,1./convergence_time_GNGD,'Linewidth',1)
xlabel('\rho value','Fontsize',11); ylabel('Convergence speed','Fontsize',11);
title('Effect of $\rho$ on algorithms','interpreter','Latex','Fontsize',11)
grid on; grid minor
legend('Benveniste','GNGD')

%% Checking performance with changing parameters - epsilon_0
epsilon_0=(0.01:0.01:2);
w0=0.9*ones(1,N); %Theoretical value
mu_0=0.1; %Initial step-size for GASS
mu=0.1; %Initial step-size for GNGD
rho=0.005;
for epsilon_index =1:length(epsilon_0)
    epsilon_index
    for iter=1:100
        
        % Generating one signal only
        h = random('Normal', 0, noise_var, N, 1);
        %Constructing MA model - Original sequence x(n)
        % We want x(n)=a1*h(n-1)+h(n)
        x=filter(b,a,h); %Output of filter
        
        % 1. Benveniste algorithm
        [~,~,er_benveniste,weights_benveniste] = benveniste_Method_1(h,x,mu_0,rho,num_coef);
        error_benveniste(:,iter)=er_benveniste;
        wmtx_benveniste(iter,:)=weights_benveniste(1,:);
        
        % 2. GNGD algorithm
        [~,er_GNGD,weights_GNGD] = GNGD_Method_1(h,x,mu,epsilon_0(epsilon_index),rho,num_coef);
        error_GNGD(:,iter)=er_GNGD;
        wmtx_GNGD(iter,:)=weights_GNGD(1,:);
    end
    w0_benveniste = mean(wmtx_benveniste,1);
    w0_GNGD = mean(wmtx_GNGD,1);
    
    w_error_benveniste = abs(w0_benveniste - w0);
    w_error_GNGD = abs(w0_GNGD - w0);
    c_s_ben = find( w_error_benveniste<threshold);
    if isempty(c_s_ben)
        c_s_ben=N;
    end
    convergence_time_ben(epsilon_index)= min (c_s_ben);
    c_s_GNGD = find( w_error_GNGD<threshold);
    if isempty(c_s_GNGD)
        c_s_GNGD=N;
    end
    convergence_time_GNGD (epsilon_index)= min (c_s_GNGD);
end

%%
figure(6);
subplot(1,3,2)
plot(epsilon_0,1./convergence_time_GNGD,'Linewidth',1)
xlabel('$\epsilon_{0}$ value','interpreter','Latex','Fontsize',11); ylabel('Convergence speed','Fontsize',11);
title('Effect of $\epsilon_{0}$ on GNGD algorithm','interpreter','Latex','Fontsize',11)
grid on; grid minor
ylim([0 0.01])
%% Checking performance with changing parameters - mu
epsilon_0=0.01;
w0=0.9*ones(1,N); %Theoretical value
mu_0=0.1; %Initial step-size for GASS
mu=(0.01:0.01:5); %Initial step-size for GNGD
size(mu)
rho=0.005;
convergence_time_ben=zeros(1,length(mu));
convergence_time_GNGD=zeros(1,length(mu));
for mu_index =1:length(mu)
    mu_index
for iter=1:100
    % Generating one signal only
    h = random('Normal', 0, noise_var, N, 1);
    %Constructing MA model - Original sequence x(n)
    % We want x(n)=a1*h(n-1)+h(n)
    x=filter(b,a,h); %Output of filter
       
    % 1. Benveniste algorithm
    [~,~,er_benveniste,weights_benveniste] = benveniste_Method_1(h,x,mu(mu_index),rho,num_coef);
    error_benveniste(:,iter)=er_benveniste;
    wmtx_benveniste(iter,:)=weights_benveniste(1,:);
        
    % 2. GNGD algorithm
    [~,er_GNGD,weights_GNGD] = GNGD_Method_1(h,x,mu(mu_index),epsilon_0,rho,num_coef);
    error_GNGD(:,iter)=er_GNGD;
    wmtx_GNGD(iter,:)=weights_GNGD(1,:); 
end
w0_benveniste = mean(wmtx_benveniste,1);
w0_GNGD = mean(wmtx_GNGD,1);

w_error_benveniste = abs(w0_benveniste - w0);
w_error_GNGD = abs(w0_GNGD - w0);
c_s_ben = find( w_error_benveniste<threshold);
if isempty(c_s_ben)
    c_s_ben=N;
end
c_s_GNGD = find( w_error_GNGD<threshold);
if isempty(c_s_GNGD)
    c_s_GNGD=N;
end
convergence_time_ben (mu_index)= min (c_s_ben);
convergence_time_GNGD (mu_index)= min (c_s_GNGD);
end
%%
figure(6);
subplot(1,3,3)
plot(mu,1./convergence_time_ben,'Linewidth',0.8)
hold on
plot(mu,1./convergence_time_GNGD,'Linewidth',1)
legend('Benveniste','GNGD')
xlabel('$\mu$ value','interpreter','Latex','Fontsize',11); ylabel('Convergence speed','Fontsize',11);
title('Effect of $\mu$ on algorithms','interpreter','Latex','Fontsize',11)
grid on; grid minor
xlim([mu(1) mu(end)])