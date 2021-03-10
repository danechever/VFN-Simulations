function n = getRefractiveIndex(material,lambda)
%   n = getRefractiveIndex(material,lambda)
%   Returns refractive index of 'material' at wavelength lambda.
%
%   Inputs:
%       material - String indicating which material (e.g. 'CaF2')
%       lambda - Wavelength in microns
%       
%   Outputs:
%       n - (Real) refractive index

    if(material == 'NBK7')
        %N-BK7 Sellmeier coefficients from Thorlabs
        %Originally from
        %citation
        B0 = 0;
        C1 = 0.00600069867;
        B1 = 1.03961212;
        C2 = 0.0200179144;
        B2 = 0.231792344;
        C3 = 103.560653;
        B3 = 1.01046945;
        
    elseif(material == 'CaF2')
        
        B0 = 0.33973;
        C1 = 0.09374;
        B1 = 0.69913;
        C2 = 21.18;
        B2 = 0.11994;
        C3 = 38.46^2;
        B3 = 4.35181;
        
    else
        error([material,'not yet implemented.']);
    end
    
    n = sqrt(1+B0+B1*lambda.^2./(lambda.^2-C1) + ...
                    B2*lambda.^2./(lambda.^2-C2) + ...
                    B3.*lambda.^2./(lambda.^2-C3));
        
