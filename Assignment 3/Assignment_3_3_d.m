%% Assignment 3.3 d)
close all; clear all ; clc

%% Loading EEG files
load EEG_Data_Assignment1.mat
Tsample = 1/fs; %Tiem resolution (seconds)
N=1200; %time samples to be used
time_ax = (0:N-1)*Tsample;

%Points to be used for time-frequency plot (Resolution - order of DFT-CLMS filter)
N_points = N;
k=0:N_points-1;

%Choosing segment of POz data
start_ind = 1000;
y = POz(start_ind:start_ind+N-1);

%Preprocessing data
%Remove mean value and trend
y = detrend(y - mean(y));

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
surf(time_ax,freq_axis,H,'LineStyle', 'none')
view(2);
c = colorbar;
c.Label.String = 'Magnitude';
c.Label.FontSize = 11;
ylim([0 100]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Time Frequency spectrum using DFT-CLMS (\gamma = 0)')

%% Implementing the DFT-CLMS algorithm with leakage
gamma=0.01;
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
ylim([0 100]);
xlabel('Sample Index'); ylabel('Frequency (Hz)');
title('Time Frequency spectrum using DFT-CLMS (\gamma = 0.01)')

figure;
plot(freq_axis,H(:,end))
xlim([0 600]);