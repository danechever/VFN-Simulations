function [phz,ptype] = vortexSelection(vtype, inpar)

    if strcmpi(vtype, 'spp')
        %***Scalar Phase Plate (Circular Baseline)***
        phz = angle(generateVortexMask(inpar.charge,inpar.coordsPP,[0 0]));
%         phz = angle(phz(:,:,ceil(inpar.numWavelengths/2)));
        ptype = 'Scalar Spiral Phase Plate';
        
    elseif strcmpi(vtype, 'spm')
        phz = generateVortexMaskKeckPrimary(inpar);%angle(makeKeckPupilInputs( inputs, initial));
        ptype = 'Segmented Primary Mirror';
        
    elseif strcmpi(vtype, 'dm')
        phz = generateDMVortex(inpar.dmBasisPath);
        ptype = 'Deformable Mirror';
    end
    
end