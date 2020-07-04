%% Assignment 3.2 a)
close all; clear all ; clc

%% General Parameters
% Time axis (n)
N=1500; %time samples to be used
n=1:N;
% Sampling frequency
fs = 1000; %Hz

%% Signal generation
%Generating circular white noise
eta = sqrt(0.05).*randn(1,N) + 1j*sqrt(0.05).*randn(1,N);

%Generating phase and frequency 
[phase,freq] = phase_generation();
phase_wrapped = wrapToPi(phase);

%Visualising frequency and phase
figure;
subplot(1,2,1)
plot(phase_wrapped)
set(gca,'YTick',-pi:pi/4:pi) 
set(gca,'YTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('Sample Index','FontSize',11)
ylabel('Phase (rad)','FontSize',11)
title('Phase \phi(n) of FM Signal','FontSize',11)
grid on
grid minor
subplot(1,2,2)
plot(freq,'Linewidth',1)
xlabel('Sample Index','FontSize',11)
ylabel('Frequency (Hz)','FontSize',11)
title('Frequency f(n) of FM Signal','FontSize',11)
grid on
grid minor

%Generating frequency-modulated (FM) signal
y = exp(1j*((2*pi*phase)/fs)) + eta;

% Finding the Circularity coefficient of data
pseudocovariance= mean(y .^ 2);
covariance = mean(abs(y) .^ 2);
circularity_coef = abs(pseudocovariance)/covariance;

figure;
scatter(real(y),imag(y))

%% Yule-Walker model - AR(1)
%AR(1) model
order=1;
[ar_coeff,~] = aryule(y,order);

%Obtaining power spectrum estimate
[h, f] = freqz(1, ar_coeff, N, fs); %h is not in dB.

figure;
subplot(1,2,1)
plot(f,pow2db(abs(h).^2),'Linewidth',1)
xlabel('Frequency (Hz)','FontSize',11)
ylabel('Magnitude (dB)','FontSize',11)
title('Power Spectrum of FM signal - AR(1) model','FontSize',11)
grid on
grid minor
xlim([0 fs/2])

%Different model orders
subplot(1,2,2)
orders=[1 2 5 10]';
for iOrder=[1 2 5 10]
    [ar_coeff,~] = aryule(y,iOrder);
    [h, f] = freqz(1, ar_coeff, N, fs);
    plot(f,pow2db(abs(h).^2),'Linewidth',1)
    hold on
end
xlabel('Frequency (Hz)','FontSize',11)
ylabel('Magnitude (dB)','FontSize',11)
title('Power Spectrum of FM signal - Different model orders','FontSize',11)
grid on
grid minor
xlim([0 fs/2])
legend(strcat('AR(',num2str(orders),')'))

%% Creating three batches of data - for the three parts of frequency
%AR(1) model for each batch
[ar_coeff_1,~] = aryule(y(1:500),1);
[ar_coeff_1_2,~] = aryule(y(1:500),5);
[ar_coeff_2,~] = aryule(y(501:1000),1);
[ar_coeff_2_2,~] = aryule(y(501:1000),5);
[ar_coeff_3,~] = aryule(y(1001:end),1);
[ar_coeff_3_2,~] = aryule(y(1001:end),5);

%Obtaining power spectrum estimates for each batch
[h_1, f_1] = freqz(1, ar_coeff_1, N, fs);
[h_1_2, f_1_2] = freqz(1, ar_coeff_1_2, N, fs);
[h_2, f_2] = freqz(1, ar_coeff_2, N, fs);
[h_2_2, f_2_2] = freqz(1, ar_coeff_2_2, N, fs);
[h_3, f_3] = freqz(1, ar_coeff_3, N, fs);
[h_3_2, f_3_2] = freqz(1, ar_coeff_3_2, N, fs);

figure;
subplot(1,3,1)
plot(f_1,pow2db(abs(h_1).^2),'Linewidth',1)
hold on
plot(f_1_2,pow2db(abs(h_1_2).^2),'Linewidth',1)
xlabel('Frequency (Hz)','FontSize',11)
ylabel('Magnitude (dB)','FontSize',11)
title('Power Spectrum of constant frequency segment','FontSize',11)
grid on
grid minor
xlim([0 fs/2])
legend('AR(1)','AR(5)')

subplot(1,3,2)
plot(f_2,pow2db(abs(h_2).^2),'Linewidth',1)
hold on
plot(f_2_2,pow2db(abs(h_2_2).^2),'Linewidth',1)
xlabel('Frequency (Hz)','FontSize',11)
ylabel('Magnitude (dB)','FontSize',11)
title('Power Spectrum of linear frequency segment','FontSize',11)
grid on
grid minor
xlim([0 fs/2])
legend('AR(1)','AR(5)')

subplot(1,3,3)
plot(f_3,pow2db(abs(h_3).^2),'Linewidth',1)
hold on
plot(f_3_2,pow2db(abs(h_3_2).^2),'Linewidth',1)
xlabel('Frequency (Hz)','FontSize',11)
ylabel('Magnitude (dB)','FontSize',11)
title('Power Spectrum of quadratic frequency segment','FontSize',11)
grid on
grid minor
xlim([0 fs/2])
legend('AR(1)','AR(5)')

%% Applying Yule-Walker every 3 time steps
%AR(1) model
order=1;
for i=1:N-2
    [ar_coeff,~] = aryule(y(i:i+2),order);
    ar_coeff_cont(:,i)=ar_coeff;
    [h, f] = freqz(1, ar_coeff, N, fs); %h is not in dB.
    h_cont(:,i)=h;
    [M I]= max(abs(h.^2));
    f_max(i) = f(I);
end
figure;
n=1:N-2;
plot(n,f_max,'Linewidth',1)
hold on
plot(freq,'Linewidth',1.5)
xlabel('Sample Index','FontSize',11)
ylabel('Frequency (Hz)','FontSize',11)
title('Frequency estimation using segment length=3','FontSize',11)
legend('AR(1) estimate','f_{0} nominal')
grid on
grid minor