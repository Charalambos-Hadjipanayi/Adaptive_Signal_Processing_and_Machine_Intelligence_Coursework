close all
clear all
clc

%% Task 1.4 b)

% Initialisation
N=1000; %Sample size
fs=1;
x = zeros(1,N);
coef = [2.76,-3.81,2.65,-0.92];
rng default %Set the random number generator to the default settings for reproducible results

for i=5:N
   x(i)=coef(1)*x(i-1) + coef(2)*x(i-2) + coef(3)*x(i-3) + coef(4)*x(i-4) + randn;    
end

%Removing first 500 samples
x=x(501:end);
p=5;
[arcoeffs e] = aryule(x,p); %p=model order and e is estimated variance of WGN input
%Note that these coefficients are inverted compared to the true
%coefficients.
%%Also these are normalised coefficients.


%% PSD estimation based on Yule-walker equations
N=length(x);

%True PSD
true_coef = [1 -coef];
[TF_true,f] = freqz(1,true_coef,[],fs);

figure;

%Generate realisation of AR process
nfft=length(TF_true);
order=(2:14);
ind=1; %Index for plotting
for i=1:length(order)
    [pxx,f] = pyulear(x,order(i),nfft*2-1,fs);
    err(i) = immse(10*log10(pxx),20*log10(abs(TF_true)));
    aic_aryule (i) = log10(err(i)) + (2*i)/N; %Akaike information criterion
    aicc_aryule (:,i) =  aic_aryule(i) + (((2*i)*(i+1))/(N-i-1)); %Corrected AIC 
    
    if order(i)==4 || order(i)==8 || order(i)==9 || order(i)==14
        figure(1);
        subplot(2,2,ind)
        plot(f,20*log10(abs(TF_true)),'Linewidth',1)
        hold on
        plot(f,10*log10(pxx),'Linewidth',1)
        grid on
        xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
        ylabel('PSD (dB)','FontSize',11)
        xlim([0 0.5])
        title(['PSD estimate by AR(' num2str(order(i)) ') model (N=500)'])
        legend(['True PSD'],['order ' num2str(order(i))])
        ind=ind+1;
    end
end

figure;
plot(order, err,'Linewidth',1)
hold on
[m index]=min(err);
stem(order(index),err(index),'r','Linewidth',1,'LineStyle','none')
xlabel('Model order','FontSize',11);
ylabel('MSE (dB)','FontSize',11)
title('MSE vs model order(N=500)','FontSize',11)
ylim([min(err) max(err)])

figure;
plot(order,aic_aryule)
hold on
plot(order,aicc_aryule)
legend('AIC criterion','AICc criterion')
xlabel('Model order')
