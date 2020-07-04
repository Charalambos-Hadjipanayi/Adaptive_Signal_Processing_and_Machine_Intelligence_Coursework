%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  CLMS Algorithm                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,weights] = CLMS_MA(x,d,mu,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% y = output signal, error = error vector, weights = matrix of weight evolution    

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N-1,1);
    weights = zeros(filter_order-1,N+1); %Stores weight time-evolution
   
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        x_n_input = x_n(1:length(x_n)-1); %This obtains only past values of x_n
        y(n) = weights(:,n)'*x_n_input;
        error(n)=d(n)-y(n); %Error calculation
        % weights update rule
        weights(:,n+1)=weights(:,n)+mu*conj(error(n))*x_n_input;
    end
    weights = weights(:,2:end); 
end