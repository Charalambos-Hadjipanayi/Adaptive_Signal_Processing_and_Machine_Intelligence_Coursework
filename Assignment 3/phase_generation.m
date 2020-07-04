function [phase,freq] = phase_generation()
    freq = zeros(1,1500);
    for n=1:1500
        if n<=500
            freq(n) = 100;
        elseif (n>500)&&(n<=1000)
            freq(n) = 100 + ((n-500)/2);
        elseif n <= 1500 
            freq(n) = 100 + (((n-1000)/25)^2);
        else
            freq(n)=0;
        end
    end
    phase = cumsum(freq);
end