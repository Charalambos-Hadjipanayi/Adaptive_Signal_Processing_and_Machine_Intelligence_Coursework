%% Assignment 3.2 a)
close all; clear all ; clc

%% General Parameters
% Time axis (n)
N=1500; %time samples to be used
n=1:N;
% Sampling frequency
fs = 1500; %Hz

%% Signal generation
%Generating circular white noise
eta = sqrt(0.05).*randn(1,N) + 1j*sqrt(0.05).*randn(1,N);

%Generating phase and frequency 
[phase,freq] = phase_generation();
phase_wrapped = wrapTo2Pi(phase);

%Generating frequency-modulated (FM) signal
y = exp(1j*((2*pi*phase)/fs)) + eta;

%% Implementing the CLMS algorithm - data is mostly circular
coef=1; %AR(1) process
N_points = 1024;
index=1; %Index for subplots
figure;
for mu=[0.005,0.01,0.4]
    
    %Balanced system
    [~,~,ar_CLMS] = CLMS_pred_AR_3_2(y,mu,coef);
    
    H = zeros(N_points,N);
    for k = 1:N
        % Compute power spectrum
        [h_CLMS ,f_CLMS]= freqz(1, [1; -conj(ar_CLMS(k)).'], N_points, fs);
        H(:, k) = abs(h_CLMS).^2; % Store it in a matrix
    end
    
    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    % Plot time-frequency diagram
    subplot(1,3,index)
    surf(n,f_CLMS,H,'LineStyle', 'none')
    view(2);
    c = colorbar;
    c.Label.String = 'PSD';
    c.Label.FontSize = 11;
    ylim([0 600]);
    xlabel('Sample Index'); ylabel('Frequency (Hz)');
    title(['Time Frequency spectrum using CLMS algorithm (\mu=',num2str(mu),')'])
    index=index+1 ;
end