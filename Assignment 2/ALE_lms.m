%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LMS Algorithm  for ALE            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [x_hat,error,weights] = ALE_lms(x,mu,order,D)
%----Inputs-----
%x = input signal, mu = step-size,
%order = system order (filter length),D = delay

%----Outputs-----
% x_hat = output signal, error = error vector, weights = matrix of weight evolution

N=length(x); %Assuming eta has same length as input_x
%Initialising variables
x_hat = zeros(N,1);
error=zeros(N,1);
weights = zeros(order,N+1); %Stores weight time-evolution

X = zeros(1,N+order+D-1);
X(order+D:N+order+D-1) = x;
x_n = zeros(order,N);

for n=1:N
    x_n(:,n) = X(n+order-1:-1:n)';
    x_hat(n) = x_n(:,n)'*weights(:,n);
    error(n)=X(n+order+D-1)-x_hat(n); %Error calculation
    % Weights update rule
    weights(:,n+1)=weights(:,n)+ mu*error(n)*x_n(:,n);
end
weights = weights(:,2:end);
end