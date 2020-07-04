%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LMS Algorithm                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,weights] = adaptive_lms(input_x,adapt_gain,filter_order)
    N=length(input_x); 
    %Initialising variables
    x_hat = zeros(N,1);
    error=zeros(N,1);
    weights = zeros(filter_order,N+1); %Stores weight time-evolution
    
    x_n = zeros(filter_order, N);
    %Obtaining x_n
    for k = 1: filter_order
        x_n(k, :) = [zeros(1, k), input_x(1: N - k)'];
    end
    
    for n=1:N
        x_hat(n) = weights(:,n)'*x_n(:,n);
        %x_hat(n) corresponds to x(n+1)
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %LMS update rule
        weights(:,n+1)=weights(:,n)+adapt_gain*error(n)*x_n(:,n);
    end
    weights = weights(:,2:end); %Discarding the first term.
end