close all
clear all
clc

%% Definition of signals
L=10000; %Signal length
noise = randn(L,1); %White Gaussian noise generated

time_lag1=(-(L-1):L-1);
index = find(time_lag1==0);
acf1 = zeros(1,(2*L-1));
acf1(index) = 1; %Variance of the noise

figure(1);
subplot(1,2,1)
plot(time_lag1,acf1,'Linewidth',1);
grid on
xlabel('Time lag (sample)','FontSize',11)
ylabel('ACF','FontSize',11)
title({'Theoretical ACF','of White Gaussian Noise'},'FontSize',11)

T=[0:L-1];
sine_wave = sin(T./100);
figure(3);
plot(sine_wave);

%%
%[acf2, time_lag2] = xcorr(sine_wave, 'unbiased'); %obtaining ACF estimate of WGN
acf2 = (1/2).*cos(time_lag1./100);
time_lag2= time_lag1;

figure(2);
subplot(1,2,1)
plot(time_lag2,acf2,'Linewidth',1);
grid on
xlabel('Time lag (sample)','FontSize',11)
ylabel('ACF','FontSize',11)
title({'Theoretical ACF','of Sinusoidal Signal'},'FontSize',11)

%% Definition 1 of PSD

%Steps: 
%1) first ifftshift to shift Inverse transform of ACF(0) to beginning of array
%2) fft brings you back to ACF
%3) then fftshift is applied to get correct PSD.
%4) get either real parts of PSD (modulus will also make them positive)

psd1 = fftshift(fft(ifftshift(acf1)));
psd2 = fftshift(fft(ifftshift(acf2)));

fs=1; %Sampling frequency for normalisation
n1=length(psd1);

%Definition of frequency axis
freqAxis1 = (-n1/2:n1/2-1)*(fs/n1);

%% Definition 2 of PSD

psd_estimate_noise1 = fftshift(fft(noise,n1));
psd_estimate_noise1=(abs(psd_estimate_noise1).^2)./L;
avg_def2 = mean(psd_estimate_noise1)

psd_estimate_2 = fftshift(fft(sine_wave,n1));
psd_estimate_2=(abs(psd_estimate_2).^2)./L;

n2=length(psd_estimate_noise1);
%Definition of frequency axis
freqAxis2=(-n2/2:n2/2-1)*(fs/n2);

%% Comparison of definitions - using plots
figure(1);
subplot(1,2,2)
plot(freqAxis2,psd_estimate_noise1,'c','Linewidth',1)
hold on
yline(mean(psd_estimate_noise1),'b','Linewidth',1.2);
plot(freqAxis1,real(psd1),'r--','Linewidth',1.5) %Take real parts since imaginary parts are very small.
grid on
title({'Periodogram of WGN using','the two definitions of PSD'},'FontSize',11)
ylabel('PSD','FontSize',11)
xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
legend('Definition 1','Mean of Definition 1','Definition 2','FontSize',9)

figure(2);
subplot(1,2,2)
plot(freqAxis1,real(psd2),'b','Linewidth',1) %Take real parts since imaginary parts are very small.
hold on
plot(freqAxis2,psd_estimate_2,'r','Linewidth',1)
grid on
title({'Periodogram of sinusoidal mixture','using the two definitions of PSD'},'FontSize',11)
ylabel('PSD','FontSize',11)
xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
legend('Definition 1','Definition 2','FontSize',11)