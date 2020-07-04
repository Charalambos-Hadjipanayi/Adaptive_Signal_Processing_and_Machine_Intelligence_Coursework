%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            CLMS Algorithm  - AR prediction         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,h] = CLMS_pred_AR(input_x,adapt_gain,filter_order)
    
    N=length(input_x); 
    x_n = zeros(filter_order, N);
    %Obtaining x_n
    for k = 1: filter_order
        x_n(k, :) = [zeros(1, k), input_x(1: N - k)];
    end
    
    %Initialising variables
    error=zeros(1,N);
    h = zeros(filter_order,N+1); %Stores weight time-evolution
    x_hat = zeros(1,N);
    
%     Setting h(0) to 1.
    h(:,1) = ones(filter_order,1);
    
    for n=1:N
        x_hat(n) = h(:,n)'*x_n(:,n);
        %CLMS error
        error(n)=input_x(n)-x_hat(n); %Error calculation
        %CLMS update rule
        h(:,n+1)=h(:,n)+adapt_gain*conj(error(n))*x_n(:,n);
    end
    h = h(:,2:end); %Discarding the first term.
end