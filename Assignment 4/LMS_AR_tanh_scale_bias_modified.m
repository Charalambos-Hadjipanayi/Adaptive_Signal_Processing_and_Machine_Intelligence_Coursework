%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    LMS Algorithm - using tanh activation function  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adding scale and bias

%Function Definition
function [x_hat,error,weights] = LMS_AR_tanh_scale_bias_modified(input_x,mu,filter_order,scale)

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
       
    for n=1:N
        s = weights(2:K,n)'*aug_x_n(2:K,n);
        x_hat(n) = weights(1,n)+scale*tanh(s); %Applying scaled tanh activation function
        %LMS error
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %LMS update rule
        weights(2:K,n+1)=weights(2:K,n)+(mu*scale*(1-(tanh(s)^2))*error(n)).*aug_x_n(2:K,n);
        weights(1,n+1)=weights(1,n)+ mu*scale*error(n)*aug_x_n(1,n);
    end
    weights = weights(:,2:end); %Discarding the first term.
end