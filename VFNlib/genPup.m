function PUPIL = genPup(ptype, inpar)
    
    if strcmpi(ptype, 'keck')
        PUPIL = makeKeckPupil(2*inpar.apRad, inpar.N );
    end
end