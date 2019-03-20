function MFD = getMFD( fiber_props, lambda )
%MFD = getMFD( fiber_props )
%   Calculates the mode field diameter based on model. 
%
%   Inputs:
%       fiber_props: structure of fiber properties 
%   Outputs:
%       MFD: Mode field diameter (meters)

    if(~isfield(fiber_props,'NA'))
        NA = sqrt(fiber_props.n_core^2-fiber_props.n_clad^2); 
    else
        NA = fiber_props.NA;
    end

    % Normalized frequency 
    V = 2*pi*fiber_props.core_rad/lambda*NA; 
    % Mode field diameter
    MFD = 2*fiber_props.core_rad*(0.65 + 1.619*V.^-1.5 + 2.879*V.^-6);


end