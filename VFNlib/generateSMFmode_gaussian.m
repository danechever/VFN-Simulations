function mode = generateSMFmode_gaussian( MFD_samples, coords )
%mode = generateSMFmode_gaussian( MFD_samples, coords )
%   Generates a Gaussian fiber mode with mode field diameter given in
%   samples. 
%
%   Inputs:
%       MFD_samples: Mode field diameter in samples 
%   Outputs:
%       mode: 2D array with normalized mode field 

    mode = sqrt(2/(pi*(MFD_samples/2)^2))* ...
        exp(-(coords.RHO/(MFD_samples/2)).^2);

end

