%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            LMS Algorithm  - AR prediction          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,weights] = LMS_pred_AR(input_x,mu,filter_order)

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
    
    %Initialising variables
    error=zeros(1,N);
    weights = zeros(filter_order,N+1); %Stores weight time-evolution
    x_hat = zeros(1,N);
       
    for n=1:N
        x_hat(n) = weights(:,n)'*x_n(:,n);
        %LMS error
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %LMS update rule
        weights(:,n+1)=weights(:,n)+mu*error(n)*x_n(:,n);
    end
    weights = weights(:,2:end); %Discarding the first term.
end