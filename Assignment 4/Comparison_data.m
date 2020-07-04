close all
clear all
clc

%% Comparing Zero-mean LMS and non-zero mean LMS

load 'w.mat';
load 'w_bias.mat';

figure;
subplot(1,2,1)
for i=1:size(weights,1)
    plot(weights(i,:),'Linewidth',1)
    hold on
end
xlabel('Sample index','FontSize',11)
ylabel('Weight amplitude','FontSize',11)
title('Weight evolution for zero mean data','FontSize',11)
ylim([-0.015 0.025])
labels_1 = (0:size(weights,1)-1);
legend(strcat('w_{',num2str(labels_1'),'}'),'Location','best')
grid on 
grid minor

subplot(1,2,2)
for i=1:size(weights_bias,1)
    plot(weights_bias(i,:),'Linewidth',1)
    hold on
end
xlabel('Sample index','FontSize',11)
ylabel('Weight amplitude','FontSize',11)
title('Weight evolution for non-zero mean data','FontSize',11)
ylim([-0.025 0.025])
labels_2 = (0:size(weights_bias,1)-1);
legend(strcat('w_{',num2str(labels_2'),'}'),'Location','best')
grid on 
grid minor

%% Comparing errors
load 'e.mat';
load 'e_bias.mat';

figure;
plot(pow2db((error.^2)./max(error.^2)))
hold on
plot(pow2db((error_bias.^2)./max(error_bias.^2)))