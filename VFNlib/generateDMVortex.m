function dmVortex = generateDMVortex(dmBasisPath)

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
dm.NdmPad = Narr;

%% Apply a vortex in voltage space 

% Calculate vortex shaped voltage map
% peakH = charge*wvl/4; % The surface height at pi for a symmetric vortex 
% [X,Y] = meshgrid(-dm.Nact/2:dm.Nact/2-1);
% dm.V = (peakH./dm.VtoH).*atan2(Y,X); 
% 
% DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);
%% Apply vortex using least squares approximation

[X2,Y2] = meshgrid(-Narr/2:Narr/2-1);
phzDesired = charge*atan2(Y2,X2); 
surfDesired = phzDesired*wvl/4/pi;
surfDesiredVec = surfDesired(:);


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

%dm = falco_update_dm_gain_map(dm);

DMsurf = falco_gen_dm_surf(dm, dm.dx, Narr);% Generate the DM surface 

E = exp(1i*4*pi*DMsurf/wvl);
E = pad_crop(E, 2048);

dmVortex = angle(E);
