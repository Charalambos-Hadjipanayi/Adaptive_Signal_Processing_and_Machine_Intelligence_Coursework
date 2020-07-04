%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LMS Algorithm - using Pre-training of coefficients %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding scale and bias + Pretraining

%Function Definition
function [x_hat,error,weights] = LMS_AR_pretraining(input_x,mu,filter_order,scale,w_init)

%----Inputs----- 
% input_x = input signal, mu = step-size, L = system order (filter length)

%----Outputs-----
% x_hat = output signal, error = error vector, weights = matrix of weight evolution  

    N=length(input_x); 
    x_n = zeros(filter_order, N);
    %Obtaining x_n
    for k = 1: filter_order
        x_n(k, :) = [zeros(1, k), input_x(1: N - k)'];
    end
    
    %Constructing augmented input - to take care of bias
    aug_x_n = [ones(1, size(x_n, 2)); x_n];
    K=size(aug_x_n,1); %This is the new filter_order
    
    %Initialising variables
    error=zeros(1,N);
    weights = zeros(K,N+1); %Stores weight time-evolution
    x_hat = zeros(1,N);
    
    %Using pre-training weights
    weights(:,1)=w_init;
       
    for n=1:N
        s = weights(:,n)'*aug_x_n(:,n);
        x_hat(n) = scale*tanh(s); %Applying scaled tanh activation function
        %LMS error
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %LMS update rule
        weights(:,n+1)=weights(:,n)+mu*scale*(1-(tanh(s).^2))*error(n)*aug_x_n(:,n);
    end
    weights = weights(:,2:end); %Discarding the first term.
end