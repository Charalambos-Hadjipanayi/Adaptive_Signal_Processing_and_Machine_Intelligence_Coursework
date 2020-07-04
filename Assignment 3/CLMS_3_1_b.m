%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  CLMS Algorithm                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,weights] = CLMS_3_1_b(x,d,mu,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% y = output signal, error = error vector, weights = matrix of weight evolution    

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N-1,1);
    weights = zeros(filter_order,N+1); %Stores weight time-evolution
    
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        y(n) = weights(:,n)'*x_n;
        error(n)=d(n)-y(n); %Error calculation
        % weights update rule
        weights(:,n+1)=weights(:,n)+mu*conj(error(n))*x_n;
    end
    weights = weights(:,2:end); 
end