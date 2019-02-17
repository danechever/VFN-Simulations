function mode = generateSMFmode_gaussian( diam_samples, coords )
%mode = generateSMFmode_gaussian( diam_samples, coords )
%   Generates a Gaussian fiber mode with mode field diameter given in
%   samples. 
%
%   Inputs:
%       diam_samples: Mode field diameter in samples 
%   Outputs:
%       mode: 2D array with normalized mode field 

    mode = sqrt(2/(pi*(diam_samples/2)^2))* ...
        exp(-(coords.RHO/(diam_samples/2)).^2);

end

