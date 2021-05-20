function GDM = buildDMbasis(dm,Narr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    disp('Building DM basis...')

    NactTotal = dm.Nact^2 ;
    dmSurfCube = zeros(Narr,Narr,NactTotal);
    for actIndex = 1:NactTotal
        dm.V = zeros(dm.Nact);
        dm.V(actIndex) = 1; 
        dmSurfCube(:,:,actIndex) = falco_gen_dm_surf(dm,dm.dx,Narr);
    end
    GDM = reshape(dmSurfCube,[Narr*Narr,NactTotal]);

    disp('Done.')
    
end

