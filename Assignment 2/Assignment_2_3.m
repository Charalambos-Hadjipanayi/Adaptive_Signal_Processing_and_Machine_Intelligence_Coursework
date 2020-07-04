%% Assignment 2.3
close all; clear all ; clc

%Number of samples to use
N=10000;

%Generating clean signal
omega_0 = 0.01*pi;
n=(0:N-1)';
x = sin(omega_0.*n);

%Filter coefficients for coloured noise
a=1;
b=[1 0 0.5];

%% ==================================================================================
%% Part a) - Verifying the minimum delta
%% ==================================================================================
mu=0.01; %LMS learning rate
filter_order = 5; %LMS filter order
N_iter=100; %Number of iterations
delay=(1:4);
figure;
MSPE_2 = zeros(length(delay),1);
figure(1);
for i=1:length(delay)
    MSPE_delay = zeros(N_iter,1);
    x_approx_1 = zeros(N_iter,N);
    subplot(1,length(delay),i)
    for iter = 1:N_iter
        v = randn(N,1); %white noise
        h=filter(b,a,v); %obtaining coloured noise
        s = x+h; %Noise-corrupted signal
        if iter==1
            plot(s,'b','Linewidth',1)
        else
            plot(s,'b','Linewidth',1,'HandleVisibility','off')
        end
        hold on
        [x_hat,~,~] = ALE_lms(s,mu,filter_order,delay(i));
        x_approx_1(iter,:)=x_hat;
        MSPE_delay(iter) = mean((x(500:end)-x_hat(500:end)).^2);
    end
    for k=1:N_iter
        if k==1
            plot(x_approx_1(k,:),'r','Linewidth',1)
        else 
            plot(x_approx_1(k,:),'r','Linewidth',1,'HandleVisibility','off')
        end
    end
    
    plot(x,'k','Linewidth',1)
    grid on
    grid minor
    ylim([-6 6])
    xlabel('Sample Index','FontSize',11); ylabel('Amplitude','FontSize',11);
    MSPE_2(i) = mean(MSPE_delay);
    title(['\Delta=',num2str(delay(i)),' and MSPE=',num2str(MSPE_2(i))])
    legend('Noise-corrpupted','ALE estimate','clean')
end


%% ==================================================================================
%% Part b) 
%% ==================================================================================

%% 1. Effects of delta on MSPE
mu=0.01; %LMS learning rate
filter_order = [5:5:20]; %LMS filter order
N_iter=100; %Number of iterations
delay = (1:25);
MSPE_1 = zeros(length(delay),1);
x_hat_plot = zeros(3,N);

index=1;
figure(2);
subplot(1,3,1)
for k=1:length(filter_order)
    k
    for i=1:length(delay)
        MSPE_delay = zeros(N_iter,1);
        x_approx_1 = zeros(N_iter,N);
        for iter = 1:N_iter
            v = randn(N,1);
            h=filter(b,a,v);
            s = x+h;
            [x_hat,~,~] = ALE_lms(s,mu,filter_order(k),delay(i));
            x_approx_1(iter,:)=x_hat;
            MSPE_delay(iter) = mean((x-x_hat).^2);
        end
        if filter_order(k)==5 && (delay(i)==3 || delay(i)==25)
            x_hat_plot(index,:)=mean(x_approx_1);
            index=index+1;
        end
        MSPE_1(i) = mean(MSPE_delay);
    end
    
    plot(MSPE_1,'Linewidth',1)
    hold on
end
xlim([1 length(delay)])
xlabel('Delay \Delta','FontSize',11); ylabel('MSPE','FontSize',11)
title('Effects of \Delta on MSPE','FontSize',11)
legend('M=5','M=10','M=15','M=20')
grid on
grid minor

figure (2)
delays_plot=[3,25];
for i=1:length(delays_plot)
    subplot(1,3,i+1)
    plot(n,x_hat_plot(i,:),'Linewidth',1)
    hold on
    plot(n,x,'k','Linewidth',1)
    xlabel('Sample Index','Fontsize',11); ylabel('Amplitude','Fontsize',11)
    title(['\Delta=', num2str(delays_plot(i))],'Fontsize',11)
    legend('$\hat{x}(n)$','$x(n)$','interpreter','latex');
    grid on
    grid minor
end


% legend()
    
    
    
%% 2. Effects of model order M on MSPE
mu=0.001; %LMS learning rate
filter_order = [1:20]; %LMS filter order
N_iter=100; %Number of iterations
delay = 3; %minimum delay used
MSPE_2 = zeros(length(filter_order),1);
N=10000;
omega_0 = 0.01*pi;
n=(0:N-1)';
x = sin(omega_0.*n);
figure;
for i=1:length(filter_order)
    i
    MSPE_model_order = zeros(N_iter,1);
    x_approx_2 = zeros(N_iter,N);
    for iter = 1:N_iter
        v = randn(N,1);
        h=filter(b,a,v);
        s = x+h;
        [x_hat,~,~] = ALE_lms(s,mu,filter_order(i),delay);
        x_approx_2(iter,:)=x_hat;
        MSPE_model_order(iter) = mean((x-x_hat).^2);
    end
    MSPE_2(i) = mean(MSPE_model_order);
%     if (filter_order(i)==1 || filter_order(i)==200)
        plot(mean(x_approx_2))
        hold on
%     end
end
plot(x)

%Obtaining optimal order
[Min_M, index] = min(MSPE_2);

figure(3);
plot(MSPE_2,'Linewidth',1)
hold on
stem(index,MSPE_2(index),'Linewidth',1,'LineStyle','none')
ylim([0.3 0.55])
xlabel('Filter order M','FontSize',11); ylabel('MSPE','FontSize',11)
title('Effects of M on MSPE','FontSize',11)
grid on
grid minor

%% ==================================================================================
%% Part c) - Implementing ANC algorithm
%% ==================================================================================

mu=0.01; %LMS learning rate
filter_order = 5; %LMS filter order
N_iter=100; %Number of iterations
x_approx_3= zeros(N_iter,N);
x_approx_4 = zeros(N_iter,N);
MSPE_3 = zeros(N_iter,1);
MSPE_4 = zeros(N_iter,1);
delay = 3; %minimum delay used
figure;
for iter = 1:N_iter
    v = randn(N,1); %white noise
    h=filter(b,a,v); %obtaining coloured noise
    s = x+h; %Noise-corrupted signal
    u=0.7*h + 0.1;%Obtaining reference signal
    
    %Computing ALE algorithm
    [x_hat,~,~] = ALE_lms(s,mu,filter_order,delay);
    x_approx_3(iter,:)=x_hat;
    MSPE_3(iter) = mean((x-x_hat).^2);
    
    %Computing ANC algorithm
    [noise_est,x_hat,~] = ANC_lms(u,s,mu,filter_order);
    x_approx_4(iter,:)=x_hat;
    MSPE_4(iter) = mean((x-x_hat).^2);
    
    subplot(2,2,1)
    if iter==1 
        plot(s,'b','Linewidth',1)
    else
        plot(s,'b','Linewidth',1,'HandleVisibility','off')
    end
    hold on
    
    subplot(2,2,2)
    if iter==1
        plot(s,'b','Linewidth',1)
    else
        plot(s,'b','Linewidth',1,'HandleVisibility','off')
    end
    hold on
end
for k=1:N_iter
    subplot(2,2,1)
    if k==1
        plot(x_approx_3(k,:),'r','Linewidth',1)
    else
        plot(x_approx_3(k,:),'r','Linewidth',1,'HandleVisibility','off')
    end
end
plot(x,'k','Linewidth',1)
grid on
grid minor
ylim([-6 6])
xlabel('Sample Index','FontSize',11); ylabel('Amplitude','FontSize',11);
legend('Noise-corrupted','ALE estimate','clean')
MSPE_ALE = mean(MSPE_3);
title(['ALE with MSE = ', num2str(MSPE_ALE),' | M=5 and \Delta=3'])

for k=1:N_iter
    subplot(2,2,2)
    if k==1
        plot(x_approx_4(k,:),'r','Linewidth',1)
    else
        plot(x_approx_4(k,:),'r','Linewidth',1,'HandleVisibility','off')
    end
end
plot(x,'k','Linewidth',1)
grid on
grid minor
ylim([-6 6])
xlabel('Sample Index','FontSize',11); ylabel('Amplitude','FontSize',11);
legend('Noise-corrupted','ANC estimate','clean')
MSPE_ANC = mean(MSPE_4);
title(['ANC with MSE = ', num2str(MSPE_ANC),' | M=5'])

%Obtaining the ensemble-averaged signal
x_approx_ALE = mean(x_approx_3);
x_approx_ANC = mean(x_approx_4);

subplot(2,2,3)
plot(x_approx_ALE,'Linewidth',1)
hold on
plot(x,'k','Linewidth',1)
grid on
grid minor
xlabel('Sample index','FontSize',11); ylabel('Amplitude','FontSize',11);
legend('ALE estimate','Clean')
title('ALE algorithm','FontSize',11)

subplot(2,2,4)
plot(x_approx_ANC,'r','Linewidth',1)
hold on
plot(x,'k','Linewidth',1)
grid on
grid minor
xlabel('Sample index','FontSize',11); ylabel('Amplitude','FontSize',11);
legend('ANC estimate','Clean')
title('ANC algorithm','FontSize',11)


