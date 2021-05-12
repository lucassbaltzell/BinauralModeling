function f = getGreenwoodCF(x)
%This function returns Greenwood center frequencies for each location x
%along the normalized cochlea
%OUTPUT
%f = center frequency
%INPUT
%x = location along normalized cochlea [0 1]

f = 165.4*(10.^(2.1*x) - 0.88);
 
end