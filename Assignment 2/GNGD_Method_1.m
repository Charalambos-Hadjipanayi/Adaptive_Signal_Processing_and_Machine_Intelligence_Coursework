%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  GNGD Algorithm                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y,error,weights] = GNGD_Method_1(x,d,mu,epsilon_0,rho,filter_order)
%----Inputs----- 
%x = input signal, d = desired signal, mu = step-size, 
%L = system order (filter length),rho = learning step-size

%----Outputs-----
% y = output signal, error = error vector, weights = matrix of weight evolution    

    N=length(d); %Assuming eta has same length as input_x
    %Initialising variables
    y = zeros(N,1);
    error=zeros(N,1);
    weights = zeros(filter_order-1,N+1); %Stores weight time-evolution
    
    X = zeros(1,N+filter_order-1);
    X(filter_order:N+filter_order-1) = x;
    epsilon = (1/mu)*ones(N+1,1);
    epsilon(1)=epsilon_0;
    x_n = zeros(filter_order,N);
    beta=1;
    
    for n=1:N
        x_n(:,n) = X(n:n+filter_order-1)';
        x_n_input(:,n) = x_n(1:length(x_n(:,n))-1,n); %This obtains only past values of x_n
        x_n_error(:,n) = x_n(end,n);
        y(n) = x_n_input(:,n)'*weights(:,n);
        error(n)=d(n)-(y(n) + x_n_error(:,n)); %Error calculation
        % GNDGD weights update rule
        weights(:,n+1)=weights(:,n)+ ((beta*error(n))/(epsilon(n)+(x_n_input(:,n)')*x_n_input(:,n)))*x_n_input(:,n);
        % Updating regularisation factor epsilon
        if n>1 %since nothing is defined for n=0
            epsilon(n+1) = epsilon(n) - (rho*mu)*((error(n)*error(n-1)*(x_n_input(:,n)')*x_n_input(:,n-1))/((epsilon(n-1)+(x_n_input(:,n-1)')*x_n_input(:,n-1))^2));
        end
    end
    weights = weights(:,2:end); 
end