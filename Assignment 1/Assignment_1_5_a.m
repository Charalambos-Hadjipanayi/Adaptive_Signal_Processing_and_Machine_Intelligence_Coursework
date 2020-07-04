%1.5 - Respiratory sinus arrhythmia from RR-Intervals
clc
clear all
close all
%% ==================================================================================
% Loading RRI data 
fs=4; %From Trial extraction procedure

load 'xRRI_trial1.mat';
load 'xRRI_trial2.mat';
load 'xRRI_trial3.mat';

%% ==================================================================================
% PSD estimate of trials

%% Trial 1

% 1) Standard Periodogram
N1=length(xRRI_trial1);
[std_psd1,fshift1] = periodogram(xRRI_trial1,rectwin(N1),N1,fs,'onesided');

figure(1);
subplot(1,3,1)
plot(fshift1,pow2db(std_psd1),'Linewidth',1)
ylim([-80 0])
grid on
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
title('Trial 1 - Standard Periodogram','Fontsize',11)

% 2) Averaged periodogram
% Constructing windows - Use of pwelch() with no overlap.

%Window size 50s.
Tsample = 1/fs;
size1 = 50/Tsample; %Convert 50 seconds into sample size

[psd_50s_1,f_50s_1] = pwelch(xRRI_trial1,rectwin(size1),0,N1,fs,'onesided');

figure(2);
subplot(1,3,1)
plot(f_50s_1,pow2db(psd_50s_1),'r','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 0])
title('Trial 1 - Averaged Periodogram (50s)','Fontsize',11)

%Window size 150s.
Tsample = 1/fs;
size2 = 150/Tsample; %Convert 150 seconds into sample size

figure(3);
subplot(1,3,1)
[psd_150s_1,f_150s_1] = pwelch(xRRI_trial1,rectwin(size2),0,N1,fs,'onesided');
plot(f_150s_1,pow2db(psd_150s_1),'k','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 0])
title('Trial 1 - Averaged Periodogram (150s)','Fontsize',11)

%% Trial 2
% 1) Standard Periodogram
N2=length(xRRI_trial2);
[std_psd2,fshift2] = periodogram(xRRI_trial2,rectwin(N2),N2,fs,'onesided');

figure(1);
subplot(1,3,2)
plot(fshift2,pow2db(std_psd2),'Linewidth',1)
grid on
ylim([-80 0])
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
title('Trial 2 - Standard Periodogram','Fontsize',11)

% 2) Averaged periodogram
% Constructing windows - Use of pwelch() with no overlap.

%Window size 50s.
Tsample = 1/fs;
size1 = 50/Tsample; %Convert 50 seconds into sample size

[psd_50s_2,f_50s_2] = pwelch(xRRI_trial2,rectwin(size1),0,N2,fs,'onesided');

figure(2);
subplot(1,3,2)
plot(f_50s_2,pow2db(psd_50s_2),'r','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 0])
title('Trial 2 - Averaged Periodogram (50s)','Fontsize',11)

%Window size 150s.
Tsample = 1/fs;
size2 = 150/Tsample; %Convert 150 seconds into sample size

[psd_150s_2,f_150s_2] = pwelch(xRRI_trial2,rectwin(size2),0,N2,fs,'onesided');
figure(3);
subplot(1,3,2)
plot(f_150s_2,pow2db(psd_150s_2),'k','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 0])
title('Trial 2 - Averaged Periodogram (150s)','Fontsize',11)

%% Trial 3
% 1) Standard Periodogram
N3=length(xRRI_trial3);
[std_psd3,fshift3] = periodogram(xRRI_trial3,rectwin(N3),N3,fs,'onesided');

figure(1);
subplot(1,3,3)
plot(fshift3,pow2db(std_psd3),'Linewidth',1)
grid on
ylim([-80 4])
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
title('Trial 3 - Standard Periodogram','Fontsize',11)

% 2) Averaged periodogram
% Constructing windows - Use of pwelch() with no overlap.

%Window size 50s.
Tsample = 1/fs;
size1 = 50/Tsample; %Convert 50 seconds into sample size

[psd_50s_3,f_50s_3] = pwelch(xRRI_trial3,rectwin(size1),0,N3,fs,'onesided');

figure(2);
subplot(1,3,3)
plot(f_50s_3,pow2db(psd_50s_3),'r','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 0])
title('Trial 3 - Averaged Periodogram (50s)','Fontsize',11)

%Window size 150s.
Tsample = 1/fs;
size2 = 150/Tsample; %Convert 150 seconds into sample size

[psd_150s_3,f_150s_3] = pwelch(xRRI_trial3,rectwin(size2),0,N3,fs,'onesided');
figure(3);
subplot(1,3,3)
plot(f_150s_3,pow2db(psd_150s_3),'k','Linewidth',1)
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
grid on
ylim([-80 4])
title('Trial 3 - Averaged Periodogram (150s)','Fontsize',11)