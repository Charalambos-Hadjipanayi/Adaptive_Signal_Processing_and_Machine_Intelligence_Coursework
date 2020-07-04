%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LMS Algorithm                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,weights] = leaky_lms(input_x,adapt_gain,leakage,filter_order)
    N=length(input_x); 
    x_hat = zeros(N,1);
    error=zeros(N,1);
    weights = zeros(filter_order,N+1);
    
    x_n = zeros(filter_order, N);
    %Obtaining x_n
    for k = 1: filter_order
        x_n(k, :) = [zeros(1, k), input_x(1: N - k)'];
    end
    
    for n=1:N
        x_hat(n) = weights(:,n)'*x_n(:,n);
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %Leaky LMS update rule
        weights(:,n+1)=(1-adapt_gain*leakage)*weights(:,n)+adapt_gain*error(n)*x_n(:,n);
    end
    weights = weights(:,2:end); %Discarding the first term.
end