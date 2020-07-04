%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LMS Algorithm                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,weights] = adaptive_lms_MA(x,d,mu,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%L = system order (filter length)

%----Outputs-----
% y = output signal, error = error vector, weights = matrix of weight evolution    

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N-1,1);
    weights = zeros(filter_order-1,N+1); %Stores weight time-evolution 
    X = zeros(1,N+filter_order-1);
    X(filter_order:N+filter_order-1) = x;

    for n=1:N
        x_n = X(n:n+filter_order-1)';
        x_n_input = x_n(1:length(x_n)-1); %This obtains only past values of x_n
        y(n) = x_n_input'*weights(:,n);
        error(n)=d(n)-y(n); %Error calculation
        % weights update rule
        weights(:,n+1)=weights(:,n)+mu*error(n)*x_n_input;
    end
    weights = weights(:,2:end); 
end