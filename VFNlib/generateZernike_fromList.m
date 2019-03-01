function phz = generateZernike_fromList( noll_indices, coeffs, PUPIL_SUPPORT, apRad, coords  )
%phz = generateZernike_fromList( noll_indices, coeffs, PUPIL_SUPPORT, apRad, coords  )
%   Generates the wavefront from list of Zernikes and their coefficients 
%
%   Inputs:
%       noll_indices: array of Noll indices 
%       coeffs: Coefficients associated with noll_indices (waves rms)
%       PUPIL_SUPPORT: Support of the entrance pupil
%       apRad: Aperture radius in units of samples
%       coords: Coordinate system structure 
%   Outputs:
%       phz: 2D array containing phase map

    phz = zeros(coords.N);
    count = 1;
    for noll_index = noll_indices
        Z = generateZernike(noll_index,apRad,coords.RHO,coords.THETA);
        Z = Z/sqrt(mean(Z(logical(PUPIL_SUPPORT)).^2)); % Re-normalize (useful when pupil is not a circle)
        phz = phz + 2*pi*coeffs(count)*Z;
        count = count + 1;
    end
    
end

