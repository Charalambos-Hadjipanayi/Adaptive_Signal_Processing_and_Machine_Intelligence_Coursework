%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Obtaining trials from ECG results          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

load 'RAW2.mat';
fs=1000;
figure;
plot(data)

%% Trial 1
data1=data(6403:250057);
[xRRI_trial1,fsRRI]=ECG_to_RRI(data1,fs);
time_ax1 = (0:1:length(xRRI_trial1)-1).*(1/fsRRI);
 
figure;
plot(time_ax1, xRRI_trial1)
xlabel('Time (s)')

xRRI_trial1 = xRRI_trial1 - mean(xRRI_trial1);
xRRI_trial1 = detrend(xRRI_trial1,10);

% Save to file
save('xRRI_trial1.mat','xRRI_trial1')

%% Trial 2

data2=[data(270990:348248); data(360665:492046)];
[xRRI_trial2,fsRRI]=ECG_to_RRI(data2,fs);
time_ax2 = (0:1:length(xRRI_trial2)-1).*(1/fsRRI);
 
figure;
plot(time_ax2, xRRI_trial2)
xlabel('Time (s)')

xRRI_trial2 = xRRI_trial2 - mean(xRRI_trial2);
xRRI_trial2 = detrend(xRRI_trial2,10);

% Save to file
save('xRRI_trial2.mat','xRRI_trial2')


%% Trial 3

data3=data(501832:751945);
[xRRI_trial3,fsRRI]=ECG_to_RRI(data3,fs);
time_ax3 = (0:1:length(xRRI_trial3)-1).*(1/fsRRI);
 
figure;
plot(time_ax3, xRRI_trial3)
xlabel('Time (s)')

xRRI_trial3 = xRRI_trial3 - mean(xRRI_trial3);
xRRI_trial3 = detrend(xRRI_trial3,10);

% Save to file
save('xRRI_trial3.mat','xRRI_trial3')


