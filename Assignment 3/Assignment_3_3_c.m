%% Assignment 3.3 c)
close all; clear all ; clc

%% General Parameters
% Time axis (n)
N=1500; %time samples to be used
time_ax=0:N-1;
% Sampling frequency
fs = 1500; %Hz
%Points to be used for time-frequency plot (Resolution - order of DFT-CLMS filter)
N_points = 1500;
k=0:N_points-1;

%% Signal generation
%Generating circular white noise
eta = sqrt(0.05).*randn(N,1) + 1j*sqrt(0.05).*randn(N,1);

%Generating phase and frequency 
[phase,freq] = phase_generation();
phase_wrapped = wrapTo2Pi(phase.');

%Generating frequency-modulated (FM) signal
y = exp(1j*((2*pi*phase.')/fs)) + eta;

%% Generating the complex phasor
x_n = zeros(N_points,N);
for i=1:N
    x_n(:,i)=(1/N_points)*exp(1j*2*(i-1)*pi*k/N_points).';
end

%% Implementing the DFT-CLMS algorithm
gamma=0;
mu=1;
freq_axis = k.*fs/N_points;

[y_hat,~,DFT_coeff] = DFT_CLMS(x_n,y,mu,gamma);

H=abs(DFT_coeff);
% Remove outliers in the matrix H
medianH = 50*median(median(H));
H(H > medianH) = medianH;
% Plot time-frequency diagram
figure; 
subplot(1,2,1)
surf(time_ax,freq_axis,H,'LineStyle', 'none')
view(2);
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 11;
ylim([0 600]);
xlabel('Sample Index'); ylabel('Frequency (Hz)');
title('Time Frequency spectrum using DFT-CLMS')

subplot(1,2,2)
plot(freq_axis,H(:,end))
xlim([0 600]);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency spectrum using DFT-CLMS at final time-step')
grid on;
grid minor;

%% Implementing the DFT-CLMS algorithm with leakage
gamma=0.05;
mu=1;
freq_axis = k.*fs/N_points;

[y_hat,~,DFT_coeff] = DFT_CLMS(x_n,y,mu,gamma);

H=abs(DFT_coeff);
% Remove outliers in the matrix H
medianH = 50*median(median(H));
H(H > medianH) = medianH;
% Plot time-frequency diagram
figure; 
surf(time_ax,freq_axis,H,'LineStyle', 'none')
view(2);
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 11;
ylim([0 600]);
xlabel('Sample Index'); ylabel('Frequency (Hz)');
title('Time Frequency spectrum using Leaky DFT-CLMS (\gamma = 0.05)')

figure;
plot(freq_axis,H(:,end))
xlim([0 600]);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency spectrum using Leaky DFT-CLMS at final time-step')
grid on;
grid minor;