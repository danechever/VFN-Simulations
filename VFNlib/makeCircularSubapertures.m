function [ PUPIL ] = makeCircularSubapertures( params, Ngrid )
%makeCircularSubapertures( params, Ngrid )
%   Makes circular or elliptical sub-apertures
%
%   Inputs: 
%       params - parameter structure length = number of subapertures 
%       params(i).x0 - x offset of ith aperture
%       params(i).y0 - y offset of ith aperture
%       params(i).radius - radius of ith aperture. Circle if length=1 and
%                          ellipse if length=2
%       Ngrid - Number of samples across array 
%
%   Outputs:
%       PUPIL - 2D array with ones in subapertures and zeros elsewhere
    PUPIL = zeros(Ngrid); 
    [X,Y] = meshgrid(-Ngrid/2:Ngrid/2-1);
    
    for apIndex = 1:length(params)
        x0 = params(apIndex).x0;
        y0 = params(apIndex).y0;
        rad = params(apIndex).radius;

        if length(rad)==1
            RHO = sqrt( ((X-x0)/rad).^2 + ((Y-y0)/rad).^2 );
        elseif length(rad)==2
            RHO = sqrt( ((X-x0)/rad(1)).^2 + ((Y-y0)/rad(2)).^2 );
        else
            error('params.radius must be length 1 or 2')
        end
        
        subAp = exp(-RHO.^2000);
        PUPIL = PUPIL + subAp;
    end
    
    PUPIL(PUPIL>1) = 1; 
    
end

