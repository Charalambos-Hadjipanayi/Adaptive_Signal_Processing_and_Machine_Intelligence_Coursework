%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  ACLMS Algorithm                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,h,g] = ACLMS_3_1_b(x,d,mu,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% y = output signal, error = error vector, h and g = matrices of weight evolution  

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N-1,1);
    h = zeros(filter_order,N+1); %Stores weight time-evolution
    g = zeros(filter_order,N+1); %Stores weight time-evolution
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        y(n) = h(:,n)'*x_n + g(:,n)'*conj(x_n);
        error(n)=d(n)-y(n); %Error calculation
        % weights update rule
        h(:,n+1)=h(:,n)+mu*conj(error(n))*x_n;
        g(:,n+1)=g(:,n)+mu*conj(error(n))*conj(x_n);
    end
    h = h(:,2:end); 
    g = g(:,2:end); 
end