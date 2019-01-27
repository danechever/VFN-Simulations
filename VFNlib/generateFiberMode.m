function mode = generateFiberMode( diam_samples, coords )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

mode = sqrt(2/(pi*(diam_samples/2)^2))* ...
    exp(-(coords.RHO/(diam_samples/2)).^2);

end

