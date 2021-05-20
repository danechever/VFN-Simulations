clear ;

%% FALCO/PROPER set up 

setPaths; 

%% User inputs 

dmConfig_BMCkiloDM; 

%%----- beam parameters 
Nbeam = 500;% Number of samples across the beam 
Narr = 512; % Number of samples across the array  
NactsPerBeam = 30; % Number of actuators per beam diameter 

% Calculated parameters from inputs 
beamDiam = NactsPerBeam*dm.dm_spacing;
dm.dx = beamDiam/Nbeam;

%% Poke some actuators 

dm.V = zeros(dm.Nact); 
dm.V(10,10) = 1; 
dm.V(5,10) = -1; 
dm.V(12,12) = 0.5; 

DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);

figure(1);
imagesc(dm.V); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM voltages')

figure(2);
imagesc(DMsurf*1e9); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM surface (nm)')
