function n = getRefractiveIndex(material,lambda)
%   n = getRefractiveIndex(material,lambda)
%   Returns refractive index of 'material' at wavelength lambda.
%
%   Inputs:
%       material - String indicating which material (e.g. 'CaF2')
%       lambda - Wavelength in METERS
%       
%   Outputs:
%       n - (Real) refractive index

    lambda = lambda * 1e6;
    
    switch lower(material)
        
        case 'nbk7' %'NBK7'
            %N-BK7 Sellmeier coefficients from Thorlabs
            % Matches this source as well: https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
            B0 = 0;
            C1 = 0.00600069867;
            B1 = 1.03961212;
            C2 = 0.0200179144;
            B2 = 0.231792344;
            C3 = 103.560653;
            B3 = 1.01046945;
            
            if (lambda < 0.3) || (lambda > 2.5)
                warning('NBK7 coefficients only valid between 0.3-2.5um; update them in code if you want more wavelengths')
            end

        case 'caf2' %'CaF2'
            % Source: https://refractiveindex.info/?shelf=main&book=CaF2&page=Li
            B0 = 0.33973;
            C1 = 0.09374^2;
            B1 = 0.69913;
            C2 = 21.18^2;
            B2 = 0.11994;
            C3 = 38.46^2;
            B3 = 4.35181;
            
            if (lambda < 0.15) || (lambda > 12)
                warning('CaF2 coefficients only valid between 0.15-12um; update them in code if you want more wavelengths')
            end

        case 'baf2' %'BaF2'
            % Source: https://refractiveindex.info/?shelf=main&book=BaF2&page=Li
            B0 = 0.33973;
            B1 = 0.81070;
            C1 = 0.10065^2;
            B2 = 0.19652;
            C2 = 29.87^2;
            B3 = 4.52469;
            C3 = 53.82^2;
            
            if (lambda < 0.15) || (lambda > 15)
                warning('BaF2 coefficients only valid between 0.15-15um; update them in code if you want more wavelengths')
            end

        case 'znse' %'ZnSe'
            % Source: https://refractiveindex.info/?shelf=main&book=ZnSe&page=Connolly
            B0 = 0;
            B1 = 4.45813734;
            C1 = 0.200859853^2;
            B2 = 0.467216334;
            C2 = 0.391371166^2;
            B3 = 2.8956629;
            C3 = 47.1362108^2;

            if (lambda < 0.54) || (lambda > 18.2)
                warning('ZnSe coefficients only valid between 0.54-18.2um; update them in code if you want more wavelengths')
            end

        case 'mgf2' %'MgF2'
            % Source: https://refractiveindex.info/?shelf=main&book=ZnSe&page=Connolly
            B0 = 0;
            B1 = 0.48755108;
            C1 = 0.04338408^2;
            B2 = 0.39875031;
            C2 = 0.09461442^2;
            B3 = 2.3120353;
            C3 = 23.793604^2;
        otherwise
            error([material,'not yet implemented.']);
    end
    
    n = sqrt(1+B0+B1*lambda.^2./(lambda.^2-C1) + ...
                    B2*lambda.^2./(lambda.^2-C2) + ...
                    B3.*lambda.^2./(lambda.^2-C3));
