function [modeC,modeS,modeParams,nModes] = generateFMFmode_bessel( NA, core_rad, lambda, dx, coords )
%mode = generateFMFmode_bessel( NA, core_rad, lambda, dx, coords )
%   Generates modes of a few mode fiber. 
%
%   Inputs:
%       NA: numerical aperture of the fiber, given by sqrt(n_core^2-n_clad^2)
%       core_rad: radius of the fiber core (meters) 
%       lambda: wavelength (meters)
%       dx: sample spacing in the image plane (meters) 
%       coords: coordinates structure 
%   Outputs:
%       modeC: 3D array with normalized mode fields (cosine terms)
%       modeS: 3D array with normalized mode fields (sine terms)
%       modeParams: Full list of mode parameters 
%       nModes: Number of modes supported 

    % Normalized frequency 
    V = 2*pi*core_rad/lambda*NA; 

    [modeParams, nModes] = getModeParams(V);

    for modeIndex = 1:length(modeParams)
        
        l = modeParams(modeIndex).l;
        m = modeParams(modeIndex).m;
        u = modeParams(modeIndex).u;
        w = modeParams(modeIndex).w;
        
        % radial mode field distributions inside and outside core
        mode_core = besselj(l,u*coords.RHO*dx/core_rad)/besselj(l,u);
        mode_clad = besselk(l,w*coords.RHO*dx/core_rad)/besselk(l,w);
        
        % Combine inner and outer fields into a single array 
        mode = mode_core;
        mode(coords.RHO*dx>core_rad) = mode_clad(coords.RHO*dx>core_rad);
        mode = (-1)^(m+1).*mode;
        
        modeC_ = mode.*cos(l*coords.THETA);
        modeS_ = mode.*sin(l*coords.THETA);

        % Return normalized mode field
        modeC(:,:,modeIndex) = modeC_/sqrt(sum(modeC_(:).^2));
        modeS(:,:,modeIndex) = modeS_/sqrt(sum(modeS_(:).^2));
        
    end
    
	modeC(isnan(modeC)) = 0;
	modeS(isnan(modeS)) = 0;  
    
end

function [modeParams, nModes] = getModeParams(V)

    modeIndex = 1; 
    l = 0; 
    while true
        
        uVec = linspace(0,V,1000);
        c = constraintDiffFunc(uVec,l,V); % solutions occur at the sign flips 
        signChanges = find(c(1:end-1)<0 & c(2:end) > 0);

        % If there are no solutions for this l then we are done
        if isempty(signChanges)
            if l==0
                nModes = modeIndex-1;
            else
                nModes = modeIndex-1 + sum(cell2mat({modeParams.l}) > 0);
            end
            break;
        end
        
        % Search for exact solutions around each sign change
        for m = 1:length(signChanges)
                  
            % Search where we know there is a change of sign
            minRange = uVec(signChanges(m) - 1);
            maxRange = uVec(signChanges(m) + 1);
            
            % Find the exact point of the change of sign
            uFit = fzero(@(u) constraintDiffFunc(u,l,V),[minRange, maxRange]);
            wFit = sqrt(V.^2-uFit^2); 
            
            if isreal(wFit) 
                % Record the solution
                modeParams(modeIndex).u = uFit;
                modeParams(modeIndex).l = l;
                modeParams(modeIndex).w = wFit;        
                modeParams(modeIndex).m = m;
                modeIndex = modeIndex + 1; 
            else
                disp(['Discarding mode ',num2str(l),num2str(m),'...']) 
            end
        
        end

        l = l + 1;
        
    end

end

function c = constraintDiffFunc(u,l,V)   

    w = sqrt(V^2 - u.^2);
    c = u.*besselj(l+1,u)./besselj(l,u) - w.*besselk(l+1,w)./besselk(l,w); 
    
end