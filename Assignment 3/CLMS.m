%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  CLMS Algorithm                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,weights] = CLMS(x,mu,filter_order)
%----Inputs----- 
%x = input signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% x_hat = output signal, error = error vector, weights = matrix of weight evolution    

    N=length(x); 
    %Initialising variables
    x_hat = zeros(N,1);
    error=zeros(N-1,1);
    weights = ones(filter_order,N+1); %Stores weight time-evolution
    
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        x_hat(n) = weights(:,n)'*x_n;
        error(n)=x(n)-x_hat(n); %Error calculation
        % weights update rule
        weights(:,n+1)=weights(:,n)+mu*conj(error(n))*x_n;
    end
    weights = weights(:,2:end); 
end