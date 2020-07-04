%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  LMS Algorithm  for ANC            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [noise_Est,x_hat,weights] = ANC_lms(u,s,mu,order)
%----Inputs-----
%s = input signal,u=reference signal (adaptive filter input), mu = step-size,
%order = system order (filter length)

%----Outputs-----
% x_hat = output signal, error = error vector, weights = matrix of weight evolution

N=length(s); %Assuming eta has same length as input_x
%Initialising variables
noise_Est = zeros(N,1);
x_hat=zeros(N,1);
weights = zeros(order,N+1); %Stores weight time-evolution

U = zeros(1,N+order-1);
U(order:N+order-1) = u;
u_n = zeros(order,N);

for n=1:N
    u_n(:,n) = U(n+order-1:-1:n)';
    noise_Est(n) = u_n(:,n)'*weights(:,n);
    x_hat(n)=s(n)- noise_Est(n); %Error calculation - error = denoised signal
    % Weights update rule
    weights(:,n+1)=weights(:,n)+ mu*x_hat(n)*u_n(:,n);
end
weights = weights(:,2:end);
end