close all
clear all
clc

%% ================= Part A ==========================

%% General parameters
L=5000; %Signal length
K=2000; %L-K will give us the number of zeros added to the end
M=500; %Number of realisations
fs=80; %Normalisation frequency
f = [0.5, 1.2, 1.5]; %Frequencies of sinusoids
psd = zeros(M,2*L-1); %Matrix storing PSD estimates
t_axis = (0:K-1)/fs;

% signal = 1.5.*sin(2*pi*f(1)*t_axis)+ 2.*sin(2*pi*f(2)*t_axis)+sin(2*pi*f(3)*t_axis);
freqAxis = (-(2*L-1)/2:(2*L-1)/2-1)*(fs/(2*L-1));

%Zero padding is also used to improve fft resolution
signal=[0.4*sin(2*pi*t_axis*f(1)) + 10*sin(2*pi*t_axis*f(2)) + 0.5*sin(2*pi*t_axis*f(3)) zeros(1, L-length(t_axis))];

%% 
figure;
subplot(1,2,1)
for i=1:M
    %Generation of sinusoidal signals corrupted by noise
    corrupted_signal = signal + randn(1,L);
    %PSD estimation - using biased ACF estimate
    [acf,~] = xcorr(corrupted_signal, 'biased'); %obtaining ACF estimate of WGN
    psd(i,:) =real(fftshift(fft(ifftshift(acf))));
    plot(freqAxis,psd(i,:),'c','Linewidth',1);
    hold on
end
xlim([0 2])
grid on
%% Finding mean value of the iterations
avg_psd = mean(psd);
plot(freqAxis,avg_psd,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD','FontSize',11)
title({'PSD estimates','(different realisations and mean)'},'FontSize',11)
%% Finding std of the iterations
std_psd = std(psd);

subplot(1,2,2)
plot(freqAxis,std_psd,'r','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('Standard deviation','FontSize',11)
title({'Standard deviation','of the PSD estimate'},'FontSize',11)
xlim([0 2])
grid on

%% ================= Part B ==========================

%Plotting same estimates in dB 
figure;
subplot(1,2,1)
for i=1:M
   psd_dB(i,:) = pow2db(psd(i,:));
   plot(freqAxis,psd_dB(i,:),'c','Linewidth',1); 
   hold on 
end
xlim([0 2])
grid on
% Finding mean value of the iterations
avg_psd_dB = mean(psd_dB);
plot(freqAxis,avg_psd_dB,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD (dB)','FontSize',11)
title({'PSD estimates in dB','(different realisations and mean)'},'FontSize',11)

%% Finding std of the iterations
std_psd_dB = std(psd_dB);

subplot(1,2,2)
plot(freqAxis,std_psd_dB,'r','Linewidth',1);
xlim([0 2])
xlabel('Frequency (Hz)','FontSize',11); ylabel('Standard deviation (dB)','FontSize',11)
title({'Standard deviation in dB','of the PSD estimate'},'FontSize',11)
grid on