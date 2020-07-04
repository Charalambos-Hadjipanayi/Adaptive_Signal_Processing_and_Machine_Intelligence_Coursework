function output = clarke_tf(input)
%Input must be 3x1 vector
%Output is 3x1 vector
%Clarke Matrix C
C=sqrt(2/3).*[sqrt(2)/2   sqrt(2)/2     sqrt(2)/2;
                    1       -0.5          -0.5;
                    0     sqrt(3)/2   -sqrt(3)/2;];
%Generating output voltage
output = C*input;
end