function mask = makeLUVOIRBPupil(Nbeam,Narr)
%mask = makeLUVOIRBPupil(Nbeam,Narr)
%   Wrapper function to use the FALCO implementation of the LUVOIR-B pupil directly  

    mp.P1.full.Nbeam = Nbeam; 
    mp.P1.compact.Nbeam = Nbeam; 
    
    mp.sbp_centers = 550e-9; % dummy value 
    mp.lambda0 = 550e-9; % dummy value  
    
    mp.P1.D = 7.989; %--meters, circumscribed. The segment size is 0.955 m, flat-to-flat, and the gaps are 6 mm. %--telescope diameter [meters]. Used only for converting milliarcseconds to lambda0/D or vice-versa.
    mp.P1.wGap = 6e-3/mp.P1.D; % Fractional width of segment gaps

    mp = falco_gen_pupil_LUVOIR_B_with_phase(mp);
    mask = mp.P1.compact.mask; 
    
    mask = pad_crop(mask,Narr); 
    
end

