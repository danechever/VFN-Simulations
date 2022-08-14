%This script simulates pupil-plane scalar VFN by generating a vortex phase
%pattern using the segments of the Keck primary. Conducts analysis of the
%resulting coupling maps in a manner similar to that shown for a circular
%pupil in Ruane et. al. 2019. This can be done by either loading a 'seed'
%array of tilt and piston from a matlab file, or using the analytical
%solution that is generated by default.
%Updated 12/2021 to use MFT propagation instead of FFT

clc;
clear; 
close all; 
% addpath(['..' filesep 'VFNlib']);
% addpath(['..' filesep '..' filesep 'falco-matlab']);

addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')));
addpath(genpath(fullfile('..','..','VFN-Lab')));

py.importlib.import_module('sellmeir1');
py.importlib.import_module('triple_prism');
load pyyo;
inpar.pyyCopy = pyyCopy;
load oldetas;
inpar.etasCopy = etasCopy;
load oldetasonax
inpar.eta_onAxCopy = eta_onAxCopy;
load py1o;
inpar.py_1Copy = py_1Copy;
load ctilt;
inpar.tilt_valsCopy = tilt_valsCopy;
sellmeir1coeffs;

% Provide path where the DM Basis file is located
dmBasisPath = '/media/Data_Drive/VFN/ScalarOnPrimaryData/';
%% Input parameters 

%-- Provide regular parameters
% Define sampling info
inpar.N = 2^10; % Size of computational grid (NxN samples) 
inpar.apRad = inpar.N/2-4; % Aperture radius in samples
inpar.apDia0 = 2 * inpar.apRad;

%-- Define wavelength info
inpar.lambda0 = 2.2e-6; %central wavelength
% fracBW = 0.2; % \Delta\lambda/\lambda
inpar.fracBW = 0.1818; %\Delta\lambda/\lambda
inpar.numWavelengths = 5; % number of discrete wavelengths 
inpar.lambdas = getWavelengthVec(inpar.lambda0, inpar.fracBW, inpar.numWavelengths);% array of wavelengths (meters)

inpar.keckD = 10.949;
% inpar.lam0OverD = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

%-- Define charge of the vortex mask at each wavelength
inpar.charge = 1*ones(1,inpar.numWavelengths); % achromatic

%-- Define wavefront error at the central wavelength
  % 1) Pist, Tip, Tilt, defoc, ob.astig, ver.astig, ver.coma, hor.coma,
  % 9) ver.tref, ob.tref, spher
inpar.nolls = 4:8;
inpar.coeffs = [0.0, 0.0, 0.0, 0.0, 0.0];  % [Waves RMS]

%-- Give offsets for the vortex mask
inpar.offsetX = 0;    % [samples in the pupil-plane]
inpar.offsetY = 0; 

%-- Parameters for SMF (Thorlabs SM600 in this case)
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436;% core index (interp @220nm)
fiber_props.n_clad = 1.4381;% cladding index (interp @220nm)
fiber_props.type = 'gaussian';

%-- Define parameters for falco MFT propagator
%<> NEED TO FIGURE THESE OUT
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 2.45e-3;    %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = getMFD(fiber_props,inpar.lambda0)/(inpar.lambda0*1.4); % focal ratio of the beam at the fiber
foc = Fnum*DPup;
fprintf('Focus in use: %f [mm]\n',foc*1e3)

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

inpar.numRings = 3;
inpar.wGap = 25.4/10916*inpar.apDia0/2;

%% Wedge & ADC initial parameters

wType = 'wedge';
ptype = 'Point Source';

inpar.p = 1.0785;%inpar.lin_offset(1);%1.1404;
inpar.wedge_mat = 'CaF2';
wedge_angle = atan(890.16*(inpar.p*inpar.lambda0)/(inpar.keckD*(getRefractiveIndex(inpar.wedge_mat ,1e6*inpar.lambda0) - 1)));

%ADC Wedge Data
wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

% inpar.wedge_ADC = 0.5 * pi/180; %Thorlabs 30 arcminute CaF2 wedge
% inpar.n1_ADC = getRefractiveIndex(wedge_ADC_1, 1e6*inpar.lambda0);
% inpar.n2_ADC = getRefractiveIndex(wedge_ADC_2, 1e6*inpar.lambda0);
% inpar.n3_ADC = getRefractiveIndex(wedge_ADC_3, 1e6*inpar.lambda0);

inpar.phi1_ADC = 7.0516 * pi/180;
inpar.phi2_ADC = 3.8050 * pi/180;
inpar.phi3_ADC = 1.1465 * pi/180;

wvs = py.numpy.array(1e6*inpar.lambdas);
inpar.n0 = 0;
inpar.n1 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));
inpar.n2 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));
inpar.n3 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));

% inpar.clocking = -88.1818;%-86.3938;%-179.976/2;
inpar.clocking = deg2rad(76.3636); %

% inpar.clocking = deg2rad(77.2727); %SPP

inpar.tilt1 = 0;

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(inpar.N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
inpar.xvalsPP = coordsPP.xvals/inpar.apRad;
inpar.yvalsPP = coordsPP.yvals/inpar.apRad;
% Used in vortex mask code
inpar.xvals = inpar.xvalsPP;
inpar.yvals = inpar.yvalsPP;
% Used in wedge/ADC code
inpar.wyvals = coordsPP.yvals;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
inpar.Nwedge = Nxi;
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
inpar.xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
inpar.yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane
lambda0Fnum_meters = inpar.lambda0*foc/DPup;

% "pixel" size in pupil and image planes
coordsPP.dx  = DPup/(2*inpar.apRad);   % meters/sample
coordsFP.dx = lambda0Fnum_meters/lambda0Fnum_samp;     % meters/sample

addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
%% Create array with pupil function

PUPIL = makeCircularPupil( inpar.apRad, inpar.N );

%-- Get norm for coupling fractions (simple sum since using MFT propagator)
totalPower0 = sum(abs(PUPIL(:)));

inpar.lambdaOverD = 8;

% Used in the wedge/ADC code
inpar.lam0OverD_meters = lambda0Fnum_meters/foc; %inpar.lambda0/inpar.keckD; %lambda0Fnum_meters; 
% inpar.lam0OverD_rad = inpar.lambda0/DPup; %lambda0Fnum_meters;

figure(); 
imagesc(inpar.xvalsPP,inpar.yvalsPP,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Define pupil field
phz = generateZernike_fromList( inpar.nolls, inpar.coeffs, PUPIL, inpar.apRad, coordsPP); 

figure();
Epup = nan(inpar.N,inpar.N,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*PUPIL;
    
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
   
end
drawnow;

%% Get PSF without vortex mask
[iPSF_BB, PSF] = getPSF_mft(Epup, inpar.lambdas, foc, coordsPP, coordsFP);

figure()
imagesc(inpar.xvalsFP,inpar.yvalsFP,iPSF_BB);
axis image; 
%axis([-3 3 -3 3]);
title('broadband PSF w/o vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Broadband PSF w/o Vortex (doesn't work for keck pupil in current implementaion)

% phz = generateZernike_fromList( inpar.nolls, inpar.coeffs, PUPIL, inpar.apRad, coordsPP);
% 
% Epup = nan(inpar.N,inpar.N,inpar.numWavelengths);
% for ch = 1:inpar.numWavelengths
%     
%     Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*PUPIL;
%     
%     subplot(1,inpar.numWavelengths,ch);
%     imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epup(:,:,ch)));
%     axis image; 
%     axis([-1 1 -1 1]);
%     title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
%     colorbar; 
%     colormap(parula(256));
% end


 %% Define pupil field

%-- Load a pregenerated seed below
%initial = load('exampleFile.mat','optimum');
%disp(optSeed);

%Scalar Phase Plate (Circular Baseline)
% phz = generateVortexMask( inpar.charge, coordsPP, [0 0] );
% central_band_index = ceil(inpar.numWavelengths/2);
% phz = angle(phz(:,:,ceil(inpar.numWavelengths/2)));
% ptype = 'Scalar Spiral Phase Plate';

% Segmented Primary Vortex Phase Pattern
% phz = generateVortexMaskKeckPrimary(inpar);%angle(makeKeckPupilInputs( inputs, initial));
% ptype = 'Segmented Primary Mirror';

% Deformable Mirror Vortex Pase Pattern
% phz = generateDMVortex(dmBasisPath);
% % phz = rot90(phz);
% ptype = 'Deformable Mirror';

figure();
imagesc(phz); 
axis image; 
colorbar;
colormap(hsv);
set(gca,'ydir','normal')
title([ptype ' Vortex Phase (radians)'])
caxis([-pi pi])

%phz(:,:,ch) = angle(makeKeckPupilPhz(inputs.apDia0, inputs.N, inputs.charge));
%phz = angle(makeKeckPupilPhase(2*apRad,N,chargeC));
%phz2 = angle(makeKeckPupilField(2*apRad,N));

figure();
Epup = nan(inpar.N,inpar.N,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    
    Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*PUPIL;
    
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epup(:,:,ch)));
    axis image; 
    axis xy;
    axis([-1 1 -1 1]);
    title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(hsv);
end
drawnow;

%% Get PSF with vortex mask
addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'propcustom']);

[iPSFv_BB, PSFv] = getPSF_mft(Epup, inpar.lambdas, foc, coordsPP, coordsFP);


figure()
imagesc(inpar.xvalsFP,inpar.yvalsFP,iPSFv_BB);
axis image; 
axis([-3 3 -3 3]);
title('Broadband PSF w/ Vortex');
colorbar;%caxis([-3 0])
colormap(parula(256));
drawnow;

%% Apply Wedge Effects

figure();
if(wType == 'ADC__')
%     inpar.lam0OverD_rad = lambda0Fnum_meters;%inpar.lambda0/inpar.keckD;
    wphz = generateWedgePlate(inpar,wedge_angle,inpar.lambdas(ceil(inpar.numWavelengths/2)),'ADC__');
    
    for ch =1:inpar.numWavelengths
                 
        phzW = phz - wphz(:,:,ch);%phz + wphz(ch);
    
        %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
        Epupw(:,:,ch) = exp(1i*wphz(:,:,ch)).*PUPIL;
        Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
    
        subplot(1,inpar.numWavelengths,ch);
        imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epupw(:,:,ch)));
        axis image; 
        axis xy;
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(hsv);
    end
    
elseif(wType == 'wedge')
    
%     inpar.lam0OverD_rad = inpar.lambda0/DPup; %inpar.keckD; %lambda0Fnum_meters; 
%     disp(pyyCopy * inpar.lam0OverD_rad);
    
%     inpar.coordsPPy = coordsPP.yvals;
        

    for ch = 1:inpar.numWavelengths
        
        inpar.pyyo = -1*pyyCopy(ch);
    
%         inpar.lam0OverD_rad = inpar.lambda0/inpar.keckD;
%         inpar.lam0OverD_rad = lambda0Fnum_meters;%inpar.lambda0/inpar.keckD;
    
        %Requisite difference in phase computed here
        wphz = generateWedgePlate(inpar, 0.0175,inpar.lambdas(ch),'wedge');s
    
        phzW = phz - wphz;
    
        %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
        Epupw(:,:,ch) = exp(1i*wphz).*PUPIL;
        Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
    
        subplot(1,inpar.numWavelengths,ch);
        imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epupw(:,:,ch)));
        axis image; 
        axis xy;
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(hsv);
    end
end

%% Get PSF after wedge/ADC applied
[PSFwBB, PSFw] = getPSF_mft(Epup_Wedge, inpar.lambdas, foc, coordsPP, coordsFP);

% figure()
% for ch = 1:inpar.numWavelengths
%     
%     subplot(1,inpar.numWavelengths,ch);
%     imagesc(inpar.xvalsPP/inpar.lambdaOverD,inpar.yvalsPP/inpar.lambdaOverD,PSFw(:,:,ch));
%     axis image;
%     axis xy;
%     axis([-2 2 -2 2]);
%     title(['PSF w/ Vortex at ' num2str(inpar.lambdas(ch)*1e9) 'nm']);
%     colorbar;%caxis([-3 0])
%     colormap(parula(256));
% end
% drawnow;

figure()
imagesc(inpar.xvalsFP,inpar.yvalsFP,PSFwBB);
axis image; 
axis([-3 3 -3 3]);
title('Broadband PSF w/ Vortex');
colorbar;
colormap(parula(256));
drawnow;
%% Generate fibermode at each lambda     
%-- Iterate through wavelengths generating modes
fibmode = nan(Nxi, Nxi, inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    % Generate the fiber mode for this wavelength with proper scaling
	fibmode(:,:,ch) = generateSMFmode_mft( fiber_props, inpar.lambdas(ch), coordsFP.dx, coordsFP);
    
end

%% Calculate Throughputs
%-- Get null depth (using overlap integral)
eta_onAx = nan(1,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    eta_onAx(ch) = (abs(sum(sum(PSFw(:,:,ch).*fibmode(:,:,ch)))).^2)/totalPower0;
end

%-- Compute 2D coupling map
eta_maps = zeros(Nxi,Nxi,inpar.numWavelengths);
for ch = 1:inpar.numWavelengths
    % Compute the monochromatic coupling map (spatial units of lambda/D)
    eta_maps(:,:,ch) = generateCouplingMap( fibmode(:,:,ch), PSFw(:,:,ch), totalPower0, 5*lambda0Fnum_samp, coordsFP);

end

disp('Key Coupling Points:')
for ch = 1:inpar.numWavelengths
    fprintf('lambda = %f nm,    on-axis coup = %e,    max = %f %%\n',inpar.lambdas(ch)*1e9, eta_onAx(ch), max(eta_maps(:,:,ch),[],'all')*100);
end

%% Display broadband coupling map if applicable
% This is what the photodiode should see when broadband light is applied

if inpar.numWavelengths > 1
    % Compute BB coupling map as average of all coupling maps
      % (This assumes the input spectrum is is flat --> ie. equal power at
      % all wavelengths)
    eta_map_BB = mean(eta_maps,3);
    
    % Plot Linear Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, eta_map_BB);
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband \eta');
    colorbar;
    colormap(gray(256));
    
    % Plot Log Scale
    figure();
    imagesc(inpar.xvalsFP, inpar.yvalsFP, log10(eta_map_BB));
    axis image;
    axis([-3 3 -3 3]);
    title('Broadband log10(\eta)');
    colorbar;
    colormap(gray(256));
    
    % Print key values
    disp('---')
    fprintf('Broadband Performance:     on-axis coup = %e,    max = %f %%\n', mean(eta_onAx), max(eta_map_BB,[],'all')*100);
end

%% Characterize Null Location

%-- find the centroid of eta_maps(:,:,ch)
%-- find min in centroid
%-- find coordinates of this min in eta_maps
%-- plot the difference in these coordinates in x and y position from the
%   origin, with lambda corresponding to ch on the bottom and offset on the
%   left
% 
%-- Finds the radius of the centroid to crop the image to depending on the maximum eta in 1 slice
%   of eta_maps(:,:,ch)
Xshift = zeros(inpar.numWavelengths,1);
Yshift = zeros(inpar.numWavelengths,1);
etaw = zeros(inpar.numWavelengths, 1);

for ch = 1:inpar.numWavelengths 
    map = eta_maps(:,:,ch); %one slice of the eta_maps cube
    map_max = max(map,[],'all'); %the maximum value in cmap
    [max_ind(ch,1),max_ind(ch,2)] = find(map_max==map,1); %linear coordinates of max value
    max_rho = sqrt(((Nxi/2+1)-max_ind(ch,1))^2 + ((Nxi/2+1)-max_ind(ch,2))^2);
    
    crp = 1.5*max_rho; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    cmap_min = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (Nxi/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
   
    Xshift(ch) = Nxi/2+1-max_ind(2); %x value is wavelength, y value is offset
    Yshift(ch) = Nxi/2+1-max_ind(1);%^
    
    etaw(ch) = cmap_min;
end

% Xshiftnw = zeros(inpar.numWavelengths,1);
% Yshiftnw = zeros(inpar.numWavelengths,1);
% etanw = zeros(inpar.numWavelengths, 1);
% for ch = 1:inpar.numWavelengths 
%     mapnw = eta_maps_NOWEDGE(:,:,ch); %one slice of the eta_maps cube
%     map_maxnw = max(mapnw,[],'all'); %the maximum value in cmap
%     [max_indnw(1),max_indnw(2)] = find(map_maxnw==mapnw,1); %linear coordinates of max value
%     max_rhonw = sqrt(((inpar.N/2+1)-max_indnw(1))^2 + ((inpar.N/2+1)-max_indnw(2))^2);
%     
%     crpnw = 2*max_rhonw; %The length of one side of the cube to crop the image to
% 
%     cmapnw = map(end/2+1-floor(crpnw/2):end/2+1+floor(crpnw/2),end/2+1-floor(crpnw/2):end/2+1+floor(crpnw/2)); %the centroid
%     cmap_minnw = min(cmapnw,[],'all'); %minimum value in the centroid
%     [min_indnw(1),min_indnw(2)] = find(cmap_minnw==cmapnw); %indices of minimum value in the centroid
%     min_indnw = min_indnw + (inpar.N/2-floor(crpnw/2)); %adjust min values to reflect position in map, not cmap
%    
%     Xshiftnw(ch) = inpar.N/2-min_indnw(2); %x value is wavelength, y value is offset
%     Yshiftnw(ch) = inpar.N/2-min_indnw(1);%^
%     
%     etanw(ch) = cmap_minnw;
% end

%change lambdaOverD here 
inpar.lambdaOverD = lambda0Fnum_samp; %inpar.N/inpar.apRad/2; % lam/D in units of samples in the image plane

%-- Null shift plots for X, Y, and the trendlines that result
figure();
subplot(2,2,1);
plot(inpar.lambdas/inpar.lambda0,Xshift/inpar.lambdaOverD, '-o', 'Color', 'r');
title('Xshift');
subplot(2,2,2);
plot(inpar.lambdas/inpar.lambda0,Yshift/inpar.lambdaOverD, '-o', 'Color', 'b');
title('Yshift');

px = polyfit(inpar.lambdas/inpar.lambda0,Xshift'/inpar.lambdaOverD,1);
pxy = polyval(px,inpar.lambdas/inpar.lambda0);
subplot(2,2,3);
plot(inpar.lambdas/inpar.lambda0,pxy,'-o','Color','m')
title('X Offset Trend')
txt = ['p value: ' num2str(px)];
text(mean(inpar.lambdas/inpar.lambda0),mean(pxy),txt);

py_1 = polyfit(inpar.lambdas/inpar.lambda0,Yshift'/inpar.lambdaOverD,1);
pyy = polyval(py_1,inpar.lambdas/inpar.lambda0);
subplot(2,2,4);
plot(inpar.lambdas/inpar.lambda0,pyy,'-o','Color','g');
title('Y Offset Trend')
txt = ['p value: ' num2str(py_1)];
text(mean(inpar.lambdas/inpar.lambda0),mean(px),txt);

%% Display coupling maps
%-- Linear scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,eta_maps(:,:,ch));
%     hold on;
%     plot((max_ind(ch,1)/inpar.lambdaOverD - (Nxi/2)/inpar.lambdaOverD),(max_ind(ch,2)/inpar.lambdaOverD - (Nxi/2)/inpar.lambdaOverD),'r+');
    axis image; 
    axis([-3 3 -3 3]);
    title(['\eta at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow

%-- Log Scale
figure();
for ch = 1:inpar.numWavelengths
    subplot(1,inpar.numWavelengths,ch);
    imagesc(inpar.xvalsFP,inpar.yvalsFP,log10(eta_maps(:,:,ch)));
    axis image; 
    axis([-3 3 -3 3]);
    title(['log10(\eta) at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(gray(256));
end
drawnow