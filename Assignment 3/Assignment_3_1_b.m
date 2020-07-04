%% Assignment 3.1 b
close all; clear all ; clc

% - Wind data for 'low', 'medium' and 'high' dynamics regions.
% - Data are recorded using the Gill Instruments WindMaster, the 2D ultrasonic anemometer
% - Wind was sampled at 32 Hz and resampled at 50Hz, and the two channels correspond to the the "north" and "east" direction
% - To make a complex-valued wind signal, combine z=v_n + j v_e, where 'v' is wind speed and 'n' and 'e' the north and east directions
% - Data length = 5000 samples
% 
% In Matlab, execute e.g.
% 
% 'load('high-wind.mat') to have the two channels as column vectors

%% Loading wind data
high_wind=load('high-wind.mat');
medium_wind=load('medium-wind.mat');
low_wind=load('low-wind.mat');

%% Forming complex-valued wind signal
high_wind_complex = high_wind.v_east + 1j*high_wind.v_north; 
medium_wind_complex = medium_wind.v_east + 1j*medium_wind.v_north;
low_wind_complex = low_wind.v_east + 1j*low_wind.v_north;

%% Cicularity plots for the three wind regimes (low,medium,high)

% Finding the Circularity coefficient for high wind
pseudocovariance_HW= mean(high_wind_complex .^ 2);
covariance_HW = mean(abs(high_wind_complex) .^ 2);
circularity_coef_HW = abs(pseudocovariance_HW)/covariance_HW;

% Finding the Circularity coefficient for medium wind
pseudocovariance_MW= mean(medium_wind_complex .^ 2);
covariance_MW = mean(abs(medium_wind_complex) .^ 2);
circularity_coef_MW = abs(pseudocovariance_MW)/covariance_MW;

% Finding the Circularity coefficient for low wind
pseudocovariance_LW= mean(low_wind_complex .^ 2);
covariance_LW = mean(abs(low_wind_complex) .^ 2);
circularity_coef_LW = abs(pseudocovariance_LW)/covariance_LW;

%Plotting the circularity plots
figure;
%1. high wind
subplot(1,3,1)
scatter(real(high_wind_complex),imag(high_wind_complex),'b','.')
hold on
scatter(mean(real(high_wind_complex)),mean(imag(high_wind_complex)),150,'g','.') %Showing the center of data
xlabel('Real part - East','Fontsize',11)
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title(['High-speed wind data with |\rho| =',num2str(round(circularity_coef_HW,4))],'Fontsize',11) 
mean_high_real = mean(real(high_wind_complex));
mean_high_im = mean(imag(high_wind_complex));
var_high_real = var(real(high_wind_complex)); 
var_high_im = var(imag(high_wind_complex));
diff_var_high = abs(var_high_real - var_high_im);

%2. medium wind
subplot(1,3,2)
scatter(real(medium_wind_complex),imag(medium_wind_complex),'r','.')
hold on
scatter(mean(real(medium_wind_complex)),mean(imag(medium_wind_complex)),150,'g','.') %Showing the center of dataxlabel('Real part - East','Fontsize',11)
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title(['Medium-speed wind data with |\rho| =',num2str(round(circularity_coef_MW,4))],'Fontsize',11)
mean_medium_real = mean(real(medium_wind_complex));
mean_medium_im = mean(imag(medium_wind_complex));
var_medium_real = var(real(medium_wind_complex)); 
var_medium_im = var(imag(medium_wind_complex));
diff_var_medium = abs(var_medium_real - var_medium_im);

%3. low wind
subplot(1,3,3)
scatter(real(low_wind_complex),imag(low_wind_complex),'k','.')
hold on
scatter(mean(real(low_wind_complex)),mean(imag(low_wind_complex)),150,'g','.') %Showing the center of data
xlabel('Real part - East','Fontsize',11)
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title(['Low-speed wind data with |\rho| =',num2str(round(circularity_coef_LW,4))],'Fontsize',11)
mean_low_real = mean(real(high_wind_complex));
mean_low_im = mean(imag(high_wind_complex));
var_low_real = var(real(low_wind_complex));
var_low_im = var(imag(low_wind_complex));
diff_var_low = abs(var_low_real - var_low_im);

%% Using CLMS and ACLMS for prediction of signal
mu=[0.001,0.005,0.1]; %Learning rate for high,medium,low
filter_length = (1:30); %Filter orders to be tested
%Initialising matrices
error_CLMS_high =zeros(length(filter_length), length(high_wind_complex));
error_CLMS_medium =zeros(length(filter_length), length(medium_wind_complex));
error_CLMS_low =zeros(length(filter_length), length(low_wind_complex));
error_ACLMS_high =zeros(length(filter_length), length(high_wind_complex));
error_ACLMS_medium =zeros(length(filter_length), length(medium_wind_complex));
error_ACLMS_low =zeros(length(filter_length), length(low_wind_complex));


x1=[0; high_wind_complex(1:end-1)];
x2=[0; medium_wind_complex(1:end-1)];
x3=[0; low_wind_complex(1:end-1)];


for order_ind=1:length(filter_length)
    order_ind
    % Implementing the CLMS algorithm
    [xhat_1,error_1,~] = CLMS_3_1_b(x1,high_wind_complex,mu(1),order_ind);
    error_CLMS_high(order_ind,:) = abs(error_1).^2;
    [xhat_2,error_2,~] = CLMS_3_1_b(x2,medium_wind_complex,mu(2),order_ind);
    error_CLMS_medium(order_ind,:) = abs(error_2).^2;
    [xhat_3,error_3,~] = CLMS_3_1_b(x3,low_wind_complex,mu(3),order_ind);
    error_CLMS_low(order_ind,:) = abs(error_3).^2;
    
    %Implementing the ACLMS algorithm
    [xhat_4,error_4,~,~] = ACLMS_3_1_b(x1,high_wind_complex,mu(1),order_ind);
    error_ACLMS_high(order_ind,:) = abs(error_4).^2;
    [xhat_5,error_5,~,~] = ACLMS_3_1_b(x2,medium_wind_complex,mu(2),order_ind);
    error_ACLMS_medium(order_ind,:) = abs(error_5).^2;
    [xhat_6,error_6,~,~] = ACLMS_3_1_b(x3,low_wind_complex,mu(3),order_ind);
    error_ACLMS_low(order_ind,:) = abs(error_6).^2;       
end

%Considering Mean Squared Errors
mse_CLMS_high = pow2db(mean(error_CLMS_high,2));
mse_CLMS_medium = pow2db(mean(error_CLMS_medium,2));
mse_CLMS_low = pow2db(mean(error_CLMS_low,2));
mse_ACLMS_high = pow2db(mean(error_ACLMS_high,2));
mse_ACLMS_medium = pow2db(mean(error_ACLMS_medium,2));
mse_ACLMS_low = pow2db(mean(error_ACLMS_low,2));

%Plotting Results
figure;
subplot(1,3,1)
plot(mse_CLMS_high,'Linewidth',1);
hold on
plot(mse_ACLMS_high,'Linewidth',1);
xlim([1 length(mse_CLMS_high)])
xlabel('Model order','Fontsize',11)
ylabel('MSPE (dB)','Fontsize',11)
grid on
grid minor
title(['High-speed wind data with \mu=',num2str(mu(1))],'Fontsize',11)
legend('CLMS','ACLMS')

subplot(1,3,2)
plot(mse_CLMS_medium,'Linewidth',1);
hold on
plot(mse_ACLMS_medium,'Linewidth',1);
xlim([1 length(mse_CLMS_medium)])
xlabel('Model order','Fontsize',11)
ylabel('MSPE (dB)','Fontsize',11)
grid on
grid minor
title(['Medium-speed wind data with \mu=',num2str(mu(2))],'Fontsize',11)
legend('CLMS','ACLMS')

subplot(1,3,3)
plot(mse_CLMS_low,'Linewidth',1);
hold on
plot(mse_ACLMS_low,'Linewidth',1);
xlim([1 length(mse_CLMS_low)])
xlabel('Model order','Fontsize',11)
ylabel('MSPE (dB)','Fontsize',11)
grid on
grid minor
title(['Lowm-speed wind data with \mu=',num2str(mu(3))],'Fontsize',11)
legend('CLMS','ACLMS')

%Difference between minimum values
[min_CLMS_high, index_1] = min(mse_CLMS_high);
[min_ACLMS_high,index_2] = min(mse_ACLMS_high);
perc_diff_high = abs(((min_CLMS_high - min_ACLMS_high)/min_CLMS_high) *100);

[min_CLMS_medium,index_3] = min(mse_CLMS_medium );
[min_ACLMS_medium,index_4]  = min(mse_ACLMS_medium );
perc_diff_medium = abs(((min_CLMS_medium  - min_ACLMS_medium )/min_CLMS_medium) *100);

[min_CLMS_low,index_5] = min(mse_CLMS_low);
[min_ACLMS_low,index_6] = min(mse_ACLMS_low);
perc_diff_low = abs(((min_CLMS_low - min_ACLMS_low)/min_CLMS_low) *100);


%% For chosen model orders with minimum MSE

M_high = 4;
M_medium = 5;
M_low = 6;

%Sm
MA_coef_num=1000; % Length of the moving average filter (NOTE: Length determines cut-off frequency, inversely proportional)
MA = ones(1,MA_coef_num)/MA_coef_num; % Impulse response of the moving average filter 

% Implementing the CLMS algorithm
[CLMS_high,~,~] = CLMS_3_1_b(x1,high_wind_complex,mu(1),M_high);
[CLMS_medium,~,~] = CLMS_3_1_b(x2,medium_wind_complex,mu(2),M_medium);
[CLMS_low,~,~] = CLMS_3_1_b(x3,low_wind_complex,mu(3),M_low);

%Implementing the ACLMS algorithm
[ACLMS_high,~,~] = ACLMS_3_1_b(x1,high_wind_complex,mu(1),M_high);
[ACLMS_medium,~,~] = ACLMS_3_1_b(x2,medium_wind_complex,mu(2),M_medium);
[ACLMS_low,~,~] = ACLMS_3_1_b(x3,low_wind_complex,mu(3),M_low);

%Plotting Learning curves
figure;
subplot(3,2,1)
scatter(real(high_wind_complex),imag(high_wind_complex),'.')
hold on
scatter(real(CLMS_high),imag(CLMS_high),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('High-speed wind data','Fontsize',11)
legend('original','CLMS')

subplot(3,2,2)
scatter(real(high_wind_complex),imag(high_wind_complex),'.')
hold on
scatter(real(ACLMS_high),imag(ACLMS_high),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('High-speed wind data','Fontsize',11)
legend('original','ACLMS')

subplot(3,2,3)
scatter(real(medium_wind_complex),imag(medium_wind_complex),'.')
hold on
scatter(real(CLMS_medium),imag(CLMS_medium),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('Medium-speed wind data','Fontsize',11)
legend('original','CLMS')

subplot(3,2,4)
scatter(real(medium_wind_complex),imag(medium_wind_complex),'.')
hold on
scatter(real(ACLMS_medium),imag(ACLMS_medium),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('Medium-speed wind data','Fontsize',11)
legend('original','ACLMS')

subplot(3,2,5)
scatter(real(low_wind_complex),imag(low_wind_complex),'.')
hold on
scatter(real(CLMS_low),imag(CLMS_low),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('Low-speed wind data','Fontsize',11)
legend('original','CLMS')

subplot(3,2,6)
scatter(real(low_wind_complex),imag(low_wind_complex),'.')
hold on
scatter(real(ACLMS_low),imag(ACLMS_low),'.')
ylabel('Imaginary part - North','Fontsize',11)
grid on
grid minor
title('Low-speed wind data','Fontsize',11)
legend('original','ACLMS')





