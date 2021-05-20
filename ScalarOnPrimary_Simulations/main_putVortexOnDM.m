clear ;
close all;
%% FALCO/PROPER set up 

setPaths; 

%% User inputs 

dmConfig_BMCkiloDM; 

%%----- beam parameters 
Nbeam = 256;% Number of samples across the beam 
Narr = 256; % Number of samples across the array  
NactsPerBeam = 30; % Number of actuators per beam diameter 

wvl = 2e-6; % Wavelength (m)
charge = 1; % Topological charge 

% Calculated parameters from inputs 
beamDiam = NactsPerBeam*dm.dm_spacing;
dm.dx = beamDiam/Nbeam;

% Provide path where the DM Basis file is located
dmBasisPath = '/media/Data_Drive/VFN/ScalarOnPrimaryData/';

%% Apply a vortex in voltage space 

% Calculate vortex shaped voltage map
peakH = charge*wvl/4; % The surface height at pi for a symmetric vortex 
[X,Y] = meshgrid(-dm.Nact/2:dm.Nact/2-1);
dm.V = (peakH./dm.VtoH).*atan2(Y,X); 

figure(1);
imagesc(dm.V); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM voltages')
drawnow; 

DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);

figure(2);
imagesc(DMsurf*1e6); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM surface (um) with vortex in voltage map')
drawnow; 

%%----- sanity check 
% E = exp(1i*4*pi*DMsurf/wvl);
% phz = angle(E);
% 
% figure(11);
% imagesc(phz); 
% axis image; 
% colorbar;
% colormap(hsv);
% set(gca,'ydir','normal')
% title('Phase (radians)')
% caxis([-pi pi])

%% Least squares fit 

%%----- Desired OPD 
[X2,Y2] = meshgrid(-Narr/2:Narr/2-1);
phzDesired = charge*atan2(Y2,X2); 
surfDesired = phzDesired*wvl/4/pi;
surfDesiredVec = surfDesired(:);

figure();
imagesc(surfDesired*1e6); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM surface (um) desired')

%%----- Build or load DM basis 
basisFileLabel = ['dmBasis_Nact',num2str(dm.Nact),'_NactBeam',num2str(NactsPerBeam),...
                '_Nbeam',num2str(Nbeam),'_Narr',num2str(Narr)]; % unique string so you can re-use the basis vectors 
% Prepend path to file
basisFileLabel = [dmBasisPath, basisFileLabel];
basisFileName = [basisFileLabel,'.fits'];

if(~isfile(basisFileName))
	GDM = buildDMbasis(dm,Narr);
	fitswrite(GDM,basisFileName);
else
    GDM = fitsread(basisFileName);
end

%%----- Get DM shape

soln = GDM\surfDesiredVec; % Calculates the voltages 
dm.V = reshape(soln,[dm.Nact,dm.Nact]);% Reshapes the voltages into 2D array 
DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);% Generate the DM surface

% for i = 1:2
%     soln = GDM\surfDesiredVec; % Calculates the voltages 
%     dm.V = reshape(soln,[dm.Nact,dm.Nact]);% Reshapes the voltages into 2D array 
%     DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);% Generate the DM surface
%     
%     if i == 1
%         dm = falco_update_dm_gain_map(dm);
%     end
% end


figure();
imagesc(DMsurf*1e6); 
axis image; 
colorbar;
colormap(gray);
set(gca,'ydir','normal')
title('DM surface (um) from least squares fit')

 %% Comparing cropped to normal
Es = exp(1i*4*pi*smallDMsurf/wvl);
sphz = angle(Es);

sphz = pad_crop(sphz,2048);

El = exp(1i*4*pi*largeDMsurf/wvl);
lphz = angle(El);



figure();
imagesc(sphz); 
axis image; 
colorbar;
colormap(hsv);
set(gca,'ydir','normal')
title('Small DMsurf array Phase (radians)')
caxis([-pi pi])

