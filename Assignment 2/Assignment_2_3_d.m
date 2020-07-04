%% Assignment 2.3 d
close all; clear all ; clc

%% Loading EEG data
load 'EEG_Data_Assignment2.mat';
data = Cz;

%% Constructing time axis
dt=1/fs; %In seconds
stop_time=length(data)*dt;
t = (0:dt:stop_time-dt)'; % seconds.

%% Constructing synthetic reference input
noise_var = 0.005;
h=sqrt(noise_var)*randn(length(data),1);
f0 = 50; % Sine wave frequency (hertz)
ref_signal = h + sin(2*pi*f0*t);

% spectrogram parameters
L=5*fs; %Length of windows
perc_over = 0.5; %percentage overlap
nOverlap = round(perc_over * L); %samples overlapping
nfft = 3*L; % nFFT points

%% Noise-corrupted EEG data
figure;
spectrogram(data, hanning(L),nOverlap , nfft, fs, 'yaxis');
ylim([0 60])
title('Spectrogram of noise-corrupted EEG data (Cz)')

%% Varying parameters
%Adaptive step-size for LMS
mu=[0.001,0.01,0.1];
%Model order
M=[2,15,30];
index=1;
for order_ind =1:length(M)
    for mu_ind=1:length(mu)
        %Computing ANC algorithm
        [noise_est,x_hat,~] = ANC_lms(ref_signal,data,mu(mu_ind),M(order_ind));
                
        % Spectrogram for ANC estimate
%         figure;
%         subplot(5,2,1)
%         spectrogram(data, rectwin(nWindows), nOverlap , nfft, fs, 'yaxis');
%         title('Original noise-corrupted signal');
%         ylim([0 60]);
        subplot(3,3,index)
        spectrogram(x_hat, hanning(L),nOverlap , nfft, fs, 'yaxis');
        title(['M = ',num2str(M(order_ind)),' and \mu = ', num2str(mu(mu_ind))]);
        ylim([0 60]);
        index=index+1;
    end
end

%% Choosing the optimal values
mu_opt = 0.001;
M_opt = 15;

%Computing ANC algorithm
[noise_est,x_hat,~] = ANC_lms(ref_signal,data,mu_opt,M_opt);

figure;
subplot(1,2,1)
spectrogram(data, hanning(L),nOverlap , nfft, fs, 'yaxis');
title('Noise-corrupted EEG signal');
ylim([0 60]);
subplot(1,2,2)
spectrogram(x_hat, hanning(L),nOverlap , nfft, fs, 'yaxis');
title(['De-noised EEG signal (M=',num2str(M_opt),',\mu=',num2str(mu_opt),')']);
ylim([0 60])

% Reducing DFT samples to 5 per Hz.
% For this case, dF is 0.0125 and we want dF=0.2 Hz.
N=length(data);
N2=N/16;
%Window size 10s.
size2 = 10/dt;%Convert 10 seconds into sample size
[psd_10s,f_10s] = pwelch(data,rectwin(size2),0,N2,fs,'onesided');
[psd_xhat,f_xhat] = pwelch(x_hat(500:end),rectwin(size2),0,N2,fs,'onesided');
psd_10s = pow2db(psd_10s); %Convert to dB
psd_xhat = pow2db(psd_xhat); %Convert to dB
figure;
subplot(1,2,1)
plot(f_10s,psd_10s,'Linewidth',1)
hold on
plot(f_xhat,psd_xhat,'Linewidth',1)
legend('Noise-corrupted','ANC de-noised')
xlabel('Frequency (Hz)','Fontsize',11);
ylabel('Amplitude (dB)','Fontsize',11);
title('Periodogram of noise-corrupted and denoised data','Fontsize',11);
grid on
grid minor
xlim([0 60])

%% Difference between the periodograms
subplot(1,2,2)
psd_diff = abs(psd_10s - psd_xhat);
plot(f_10s,psd_diff)
xlabel('Frequency (Hz)','Fontsize',11);
ylabel('Absolute Error (dB)','Fontsize',11);
title('Absolute error between noise-corrupted and de-noised data','Fontsize',11);
grid on
grid minor
xlim([0 60])


