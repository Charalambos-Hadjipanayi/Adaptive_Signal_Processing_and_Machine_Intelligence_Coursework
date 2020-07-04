%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  ACLMS Algorithm                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,h,g] = ACLMS_MA(x,d,mu,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%filter_order = system order (filter length)

%----Outputs-----
% y = output signal, error = error vector, h and g = matrices of weight evolution  

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N-1,1);
    h = zeros(filter_order-1,N+1); %Stores weight time-evolution
    g = zeros(filter_order-1,N+1); %Stores weight time-evolution
    X = zeros(N+filter_order-1,1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1);
        x_n_input = x_n(1:length(x_n)-1); %This obtains only past values of x_n
        y(n) = h(:,n)'*x_n_input + g(:,n)'*conj(x_n_input);
        error(n)=d(n)-y(n); %Error calculation
        % weights update rule
        h(:,n+1)=h(:,n)+mu*conj(error(n))*x_n_input;
        g(:,n+1)=g(:,n)+mu*conj(error(n))*conj(x_n_input);
    end
    h = h(:,2:end); 
    g = g(:,2:end); 
end