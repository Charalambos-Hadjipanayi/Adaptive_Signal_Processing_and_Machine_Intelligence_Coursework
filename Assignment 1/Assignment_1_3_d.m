close all
clear all
clc

%% Generation of complex-valued signal with complex noise
fs=1;
N=[20,30,40,80,90,110];

figure(1);
subplot(1,2,1)
for i=1:3
    n = 0:N(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

    % Periodogram - Spectral estimation
    dF = fs/N(i); %indexing in Hz
    % Considering 512 frequency bins (512 fft points per Hz)
    % Increasing nFFT basically does zero padding.
    dF_new=fs/512; %To increase bins to 512, dF must be 1/512;
    K=dF_new/dF;
    [pxx,f] = periodogram(x,rectwin(N(i)),round(N(i)/K),fs);
    plot(f.*1000,pow2db(pxx),'Linewidth',1)
    hold on
end
xlabel('Frequency (mHz)','FontSize',11)
ylabel('Power/frequency (dB/Hz)','FontSize',11)
grid on
title('PDS estimates of complex exponentials with noise','FontSize',11)
legend('N=20','N=30','N=40','FontSize',9)
xlim([0.2 0.4].*1000)


figure(1);
subplot(1,2,2)
for i=4:6
    n = 0:N(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

    % Periodogram - Spectral estimation
    dF = fs/N(i); %indexing in Hz
    % Considering 512 frequency bins (512 fft points per Hz)
    % Increasing nFFT basically does zero padding.
    dF_new=fs/512; %To increase bins to 512, dF must be 1/512;
    K=dF_new/dF;
    [pxx,f] = periodogram(x,rectwin(N(i)),round(N(i)/K),fs);
    plot(f.*1000,pow2db(pxx),'Linewidth',1)
    hold on
end
xlabel('Frequency (mHz)','FontSize',11)
ylabel('Power/frequency (dB/Hz)','FontSize',11)
grid on
title('PDS estimates of complex exponentials with noise','FontSize',11)
legend('N=80','N=90','N=110','FontSize',9)
xlim([0.2 0.4].*1000)