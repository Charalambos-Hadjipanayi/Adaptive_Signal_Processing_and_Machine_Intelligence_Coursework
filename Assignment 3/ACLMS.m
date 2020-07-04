%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  ACLMS Algorithm                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,h,g] = ACLMS(x,mu,filter_order)
%----Inputs----- 
%x = input signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% x_hat = output signal, error = error vector, h and g = matrices of weight evolution  

    N=length(x); %Assuming eta has same length as input_x
    %Initialising variables
    x_hat = zeros(N,1);
    error=zeros(N-1,1);
    h = ones(filter_order,N+1); %Stores weight time-evolution
    g = ones(filter_order,N+1); %Stores weight time-evolution
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        x_hat(n) = h(:,n)'*x_n + g(:,n)'*conj(x_n);
        error(n)=x(n)-x_hat(n); %Error calculation
        % weights update rule
        h(:,n+1)=h(:,n)+mu*conj(error(n))*x_n;
        g(:,n+1)=g(:,n)+mu*conj(error(n))*conj(x_n);
    end
    h = h(:,2:end); 
    g = g(:,2:end); 
end