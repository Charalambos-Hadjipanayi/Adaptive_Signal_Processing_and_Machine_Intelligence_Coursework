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

% AR Spectrum
p=(1:15);
for i=1:length(p)
    [a,err] = aryule(xRRI_trial1,p(i));
    errors (i) = err;
    coeff(i,:)=[a zeros(1,length(p)-i)];
end

figure;
plot(errors);

% [h,f] = freqz(___,n,fs)

%From coefficients, this looks like AR(10) process
[pxx1 f1]= pyulear(xRRI_trial1,10,[],fs);

%Standard periodogram method
N1=length(xRRI_trial1);
[std_psd1,fshift1] = periodogram(xRRI_trial1,rectwin(N1),N1,fs,'onesided');

figure(4);
subplot(1,3,1)
plot(fshift1,pow2db(std_psd1),'r','Linewidth',1)
hold on
plot(f1,pow2db(pxx1),'b','Linewidth',1.5)
ylim([-80 0])
grid on
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
legend('Periodogram','AR(10)')
title('Trial 1 - AR Spectrum')

%% Trial 2

% AR Spectrum
p=(1:15);
for i=1:length(p)
    [a,err] = aryule(xRRI_trial2,p(i));
    errors (i) = err;
    coeff2(i,:)=[a zeros(1,length(p)-i)];
end

figure;
plot(errors);

% [h,f] = freqz(___,n,fs)

%From coefficients, this looks like AR(11) process
[pxx2 f2]= pyulear(xRRI_trial2,11,[],fs);

%Standard periodogram method
N2=length(xRRI_trial2);
[std_psd2,fshift2] = periodogram(xRRI_trial2,rectwin(N2),N2,fs,'onesided');

figure (4);
subplot(1,3,2)
plot(fshift2,pow2db(std_psd2),'r','Linewidth',1)
hold on
plot(f2,pow2db(pxx2),'b','Linewidth',1.5)
xline(25/60,'k--','Linewidth',1);
ylim([-80 0])
grid on
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
legend('Periodogram','AR(11)','25bpm')
title('Trial 2 - AR Spectrum')

%% Trial 3

% AR Spectrum
p=(1:15);
for i=1:length(p)
    [a,err] = aryule(xRRI_trial3,p(i));
    errors (i) = err;
    coeff3(i,:)=[a zeros(1,length(p)-i)];
end

figure;
plot(errors);

% [h,f] = freqz(___,n,fs)

%From coefficients, this looks like AR(9) process
[pxx3 f3]= pyulear(xRRI_trial3,8,[],fs);

%Standard periodogram method
N3=length(xRRI_trial3);
[std_psd3,fshift3] = periodogram(xRRI_trial3,rectwin(N3),N3,fs,'onesided');

figure (4);
subplot(1,3,3)
plot(fshift3,pow2db(std_psd3),'r','Linewidth',1)
hold on
plot(f3,pow2db(pxx3),'b','Linewidth',1.5)
xline(7.5/60,'k--','Linewidth',1);
ylim([-80 0])
grid on
xlabel('Frequency(Hz)','Fontsize',11)
ylabel('Power/Frequency(dB/Hz)','Fontsize',11)
legend('Periodogram','AR(8)','7.5 bpm')
title('Trial 3 - AR Spectrum')
