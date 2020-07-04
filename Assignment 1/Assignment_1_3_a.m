close all
clear all
clc

%% Generation of signals
L=10000;
t_axis = (0:L-1);
fs=1; %Sampling frequency for normalisation

%WGN
noise =(rand(L,1)-0.5)*2; %White Gaussian noise generated

figure;
plot(t_axis,noise)
xlabel('Time (sample index)')
ylabel('Amplitude')
title('WGN vs time')

%Noisy sinusoidal signal
f=[0.1,0.4];
noisy_sine_wave = sin((2*pi*f(1)).*t_axis)+ sin((2*pi*f(2)).*t_axis) + (randn(1,L)-0.5).*2;
figure;
plot(t_axis,noisy_sine_wave)
xlabel('Time (sample index)')
ylabel('Amplitude')
title('Noisy sinusoidal signal vs time')

%Filtered WGN - by MA filter
x=(rand(L,1)-0.5)*2;
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
filtered_noise= filter(b,a,x);
figure;
plot(t_axis,filtered_noise)
xlabel('Time (sample index)')
ylabel('Amplitude')
title('Filtered WGN vs time')

%% Applying the biased and unbiased ACF estimates

%WGN
[acf_WGN_biased, timelag_WGN] = xcorr(noise, 'biased'); %obtaining ACF estimate of WGN
[acf_WGN_unbiased,~] = xcorr(noise, 'unbiased');

figure(4);
subplot(3,2,1)
plot(timelag_WGN,acf_WGN_unbiased,'Linewidth',1);
hold on
plot(timelag_WGN,acf_WGN_biased,'Linewidth',1);
grid on
xlabel('Time lag','FontSize',11)
ylabel('ACF magnitude','FontSize',11)
legend('Unbiased','Biased','FontSize',9)
title('ACF for WGN','FontSize',11)

%Noisy sinusoidal signal
[acf_sine_biased, timelag_sine] = xcorr(noisy_sine_wave, 'biased'); %obtaining ACF estimate of WGN
[acf_sine_unbiased,~] = xcorr(noisy_sine_wave, 'unbiased');

figure(4);
subplot(3,2,3)
plot(timelag_sine,acf_sine_unbiased,'Linewidth',1);
hold on
plot(timelag_sine,acf_sine_biased,'Linewidth',1);
grid on
xlabel('Time lag','FontSize',11)
ylabel('ACF magnitude','FontSize',11)
legend('Unbiased','Biased','FontSize',9)
title('ACF for noisy sinusoidal signal','FontSize',11)

%Filtered WGN
[acf_filtered_biased, timelag_filtered] = xcorr(filtered_noise, 'biased'); %obtaining ACF estimate of WGN
[acf_filtered_unbiased,~] = xcorr(filtered_noise, 'unbiased');

figure(4);
subplot(3,2,5)
plot(timelag_filtered,acf_filtered_unbiased,'Linewidth',1);
hold on
plot(timelag_filtered,acf_filtered_biased,'Linewidth',1);
grid on
xlabel('Time lag','FontSize',11)
ylabel('ACF magnitude','FontSize',11)
legend('Unbiased','Biased','FontSize',9)
title('ACF for filtered WGN','FontSize',11)

%% Obtaining the correlogram spectral estimator
%Steps: 
%1) first ifftshift to shift Inverse transform of ACF(0) to beginning of array
%2) fft brings you back to ACF
%3) then fftshift is applied to get correct PSD.
%4) get either real parts of PSD (modulus will also make them positive)

%For WGN
psd_wgn_unbiased = fftshift(fft(ifftshift(acf_WGN_unbiased)));
psd_wgn_biased =fftshift(fft(ifftshift(acf_WGN_biased)));

n=length(psd_wgn_unbiased);
%Definition of frequency axis
freqAxis = (-n/2:n/2-1)*(fs/n);

figure(4);
subplot(3,2,2)
plot(freqAxis,real(psd_wgn_unbiased),'Linewidth',1)
hold on
plot(freqAxis,real(psd_wgn_biased),'Linewidth',1)
grid on
ylabel('PSD','FontSize',11)
xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
legend('Unbiased','Biased','FontSize',9)
title('PSD for filtered WGN','FontSize',11)

%For Noisy sinusoidal signal
psd_sine_unbiased = fftshift(fft(ifftshift(acf_sine_unbiased)));
psd_sine_biased =fftshift(fft(ifftshift(acf_sine_biased)));

n=length(psd_sine_unbiased);
%Definition of frequency axis
freqAxis = (-n/2:n/2-1)*(fs/n);

figure(4);
subplot(3,2,4)
plot(freqAxis,real(psd_sine_unbiased),'Linewidth',1)
hold on
plot(freqAxis,real(psd_sine_biased),'Linewidth',1)
grid on
ylabel('PSD','FontSize',11)
xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
legend('Unbiased','Biased','FontSize',9)
title('PSD for Noisy sinusoidal signal','FontSize',11)


%For filtered WGN
psd_filtered_unbiased = fftshift(fft(ifftshift(acf_filtered_unbiased)));
psd_filtered_biased =fftshift(fft(ifftshift(acf_filtered_biased)));

n=length(psd_filtered_unbiased);
%Definition of frequency axis
freqAxis = (-n/2:n/2-1)*(fs/n);

figure(4);
subplot(3,2,6)
plot(freqAxis,real(psd_filtered_unbiased),'Linewidth',1)
hold on
plot(freqAxis,real(psd_filtered_biased),'Linewidth',1)
grid on
ylabel('PSD','FontSize',11)
xlabel('Normalised frequency (\pi rad/sample)','FontSize',11);
legend('Unbiased','Biased','FontSize',9)
title('PSD for filtered WGN','FontSize',11)
