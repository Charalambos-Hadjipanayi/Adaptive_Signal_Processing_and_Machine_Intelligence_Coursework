%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Benveniste Algorithm                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,mu,error,weights] = benveniste_Method_1(x,d,mu_0, rho,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu_0 = initial step-size, 
%rho = learning rate step-size, L = system order (filter length)

%----Outputs-----
% y = output signal, mu = vector of step-size evolution
% error = error vector, weights = matrix of weight evolution    

    N=length(d); %Assuming d has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N,1);
    weights = zeros(filter_order-1,N+1); %Stores weight time-evolution
    mu = zeros(N+1,1);
    psi = zeros(filter_order-1,N+1);
    
    X = zeros(1,N+filter_order-1);
    X(filter_order:N+filter_order-1) = x;
    
    mu(1)=mu_0;
    for n=1:N
        x_n = X(n:n+filter_order-1)';
        x_n_input = x_n(1:length(x_n)-1); %This obtains only past values of x_n
        x_n_error = x_n(end);
        y(n) = x_n_input'*weights(:,n); 
        error(n)=d(n)-(y(n)+x_n_error); %Error calculation
        % weights update rule
        weights(:,n+1)=weights(:,n)+mu(n)*error(n)*x_n_input;
        %Updating step-size mu
        mu(n+1)= mu(n) + rho * error(n) * x_n_input'*psi(:,n);
        %Computing psi term for next loop
        psi(:,n+1) = (eye(filter_order-1,filter_order-1) - mu(n) * x_n_input*(x_n_input')) * psi(:,n) + error(n) * x_n_input;
    end
    weights = weights(:,2:end); 
end