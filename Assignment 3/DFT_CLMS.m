%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  DFT-CLMS Algorithm                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function Definition
function [y_hat,error,DFT_coeff] = DFT_CLMS(x,d,mu,gamma)
%----Inputs----- 
%x=input complex phasor, d = desired signal, mu = step-size, gamma=leakage factor

%----Outputs-----
% y_hat = output signal, error = error vector, DFT_coeff = matrix of weight evolution    

    N=length(d); %Assuming eta has same length as input_x
    K=size(x,1);
    %Initialising variables
    y_hat = zeros(N,1);
    error=zeros(N,1);
    DFT_coeff = zeros(K,N+1); %Stores weight time-evolution
    
    for n=1:N
        %Generating complex phasor (Fourier basis)
        y_hat(n) = DFT_coeff(:,n)'*x(:,n);
        error(n)=d(n)-y_hat(n); %Error calculation
        % weights update rule
        DFT_coeff(:,n+1)=(1-gamma*mu)*DFT_coeff(:,n)+mu*conj(error(n))*x(:,n);
    end
    DFT_coeff = DFT_coeff(:,2:end); 
end