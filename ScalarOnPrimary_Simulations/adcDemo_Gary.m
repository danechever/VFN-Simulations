% Demo of a simple ADC implementation

clear; close all;

% Set Paths
addpath(genpath(fullfile('..\VFN\VFN-Simulations\ScalarOnPrimary_Simulations')));
addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')));
addpath(genpath(fullfile('..','..','VFN-Lab')));
addpath(genpath(fullfile('..','..','ADC')));
addpath 'C:\Users\cyrus\OneDrive - California Institute of Technology\School Stuff\Research\VFN\ADC\RayTrace'

%% Input parameters 

sellmeir1coeffs;

%-- Provide regular parameters
% Define sampling info
N = 2^10; % Size of computational grid (NxN samples) 
apRad = N/2-4; % [samples] Aperture radius

%-- Define wavelength info
lambda0 = 2.2e-6; % [m] central wavelength
fracBW = 0.1818; % [\Delta\lambda/\lambda] fractional bandwidth
numWavelengths = 5; % number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0, fracBW, numWavelengths);% [m] array of wavelengths

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
keckD = 10.949;   %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = 4.6; % focal ratio of the beam at the final focal plane
foc = Fnum*keckD;
fprintf('Focus in use: %f [mm]\n',foc*1e3)

mag = 1;
beamD = mag*keckD;

% Since using falco MFT propogator, can choose our final image sampling
lambda0Fnum_samp = 50; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

%% Wedge & ADC initial parameters

wType = 'adc';    % options: 'wedge', 'adc'

wedge_mat = 'caf2';     % wedge material ('NBK7', 'CaF2', 'BaF2', 'ZnSe')
%n_wedge = getRefractiveIndex(wedge_mat, lambda0);   % refractive index at central wavelength
%beam_dev = 1.1*lambda0/beamD;    %[rads] Beam deviation to simulate at central wavelength
wedge_angle = deg2rad(10);%atan2(beam_dev,(n_wedge-1));   % [rads] wedge angle to use
inpar.clocking = deg2rad(34.1697);

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
xvalsPP = coordsPP.xvals/apRad;
yvalsPP = coordsPP.yvals/apRad;

%-- Coordinates in the pupil plane
Nxi = im_size*lambda0Fnum_samp;
%inpar.Nwedge = Nxi;
coordsFP = generateCoordinates( Nxi );
% Helpful for plotting: get axes coords in lam/D 
xvalsFP = coordsFP.xvals/lambda0Fnum_samp;
yvalsFP = coordsFP.yvals/lambda0Fnum_samp;

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in meters at final focal plane
lambda0Fnum_meters = lambda0*foc/keckD; % [m]

% "pixel" size in pupil and image planes
coordsPP.dx = keckD/(2*apRad);   % [meters/sample]
coordsFP.dx = lambda0Fnum_meters/lambda0Fnum_samp;     % [meters/sample]

%% Create array with pupil function

PUPIL = makeCircularPupil(apRad, N);

figure();
imagesc(xvalsPP,yvalsPP,PUPIL); 
axis image;
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Get PSF without wedge
% Replicate pupil at each wavelength for this since all have the same field
Epup = repmat(PUPIL, [1,size(lambdas)]);
[~, PSF] = getPSF_mft(Epup, lambdas, foc, coordsPP, coordsFP);

figure();
pause(0.1);
for ch = 1:numWavelengths    
    subplot(1,numWavelengths,ch);
    iPSF = abs(PSF(:,:,ch)).^2; % Compute intensity (ie. image)
    imagesc(xvalsFP,yvalsFP,iPSF);
    axis image; 
    axis xy;
    axis([-3 3 -3 3]);
    title(['PSF at ',num2str(lambdas(ch)*1e9),'nm']);
    colorbar; 
    colormap(parula(256));
    drawnow;
end
%% Generate the ADC

inpar.keckD = keckD;
inpar.lambdas = lambdas;
inpar.numWavelengths = numWavelengths;

%ADC Wedge Data
wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

% inpar.wedge_ADC = 0.5 * pi/180; %Thorlabs 30 arcminute CaF2 wedge
% inpar.n1_ADC = getRefractiveIndex(wedge_ADC_1, 1e6*inpar.lambda0);
% inpar.n2_ADC = getRefractiveIndex(wedge_ADC_2, 1e6*inpar.lambda0);
% inpar.n3_ADC = getRefractiveIndex(wedge_ADC_3, 1e6*inpar.lambda0);

inpar.phi1_ADC = -1*7.0516 * pi/180;
inpar.phi2_ADC = 3.8050 * pi/180;
inpar.phi3_ADC = 1.1465 * pi/180;

wvs = py.numpy.array(1e6*lambdas);
% inpar.n0 = 1;
inpar.n1 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));
inpar.n2 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));
inpar.n3 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));

wphz = simpleGaryADC(inpar);

figure();

% [inpar.tilt_valsCopy, inpar.clocking, inpar.I] = simple_ADCOPT_DEMO(inpar);
inpar.tilt_valsCopy = simpleGaryADC(inpar);

disp('Tilt in radians at every wavelength');
disp(inpar.tilt_valsCopy);
disp('Peak to Valley Tilt across the polychromatic band in milliarcseconds');
disp(206265*1000*(inpar.tilt_valsCopy(end) - inpar.tilt_valsCopy(1)));

inpar.magfactor = mag;

inpar.coordsPP.Y = coordsPP.Y; %Y = inputs.coordsPP.Y;
inpar.coordsPP.dx = coordsPP.dx;
    
inpar.xvalsFP = xvalsFP;
inpar.yvalsFP = yvalsFP;
inpar.xvalsPP = xvalsPP;
inpar.yvalsPP = yvalsPP;
    
disp('Peak to valley tilt across the pupil in wave fractions');

for ch = 1:inpar.numWavelengths
    Y = inpar.coordsPP.Y;
    dx = inpar.coordsPP.dx;
    
    wphz = (wphz - wphz(ceil(inpar.numWavelengths/2)));
    wphz = wphz';
    for i = 1:inpar.numWavelengths
        
        gphz(:,:,i) = 2*pi*(wphz(i)/mag/inpar.lambdas(i))*Y*dx;
        
    end
       
        phzW = gphz(:,:,ch);%phz - wphz(:,:,ch);%phz + wphz(ch);
    
        %Epup(:,:,ch) = exp(1i*phz*inpar.lambda0/inpar.lambdas(ch)).*inpar.PUPIL;
        Epupw(:,:,ch) = exp(1i*gphz(:,:,ch)).*PUPIL;
        Epup_Wedge(:,:,ch) = exp(1i*phzW).*PUPIL;
        
        disp(max(max(gphz(:,:,ch)/2/pi))-min(min(gphz(:,:,ch)/2/pi)));
        
        subplot(1,inpar.numWavelengths,ch);
        imagesc(inpar.xvalsPP,inpar.yvalsPP,angle(Epupw(:,:,ch)));
        axis image; 
        axis xy;
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(inpar.lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(hsv);
end

%% Get PSF after wedge/ADC applied
[PSFwBB, PSFw] = getPSF_mft(Epup_Wedge, inpar.lambdas, foc, coordsPP, coordsFP);

figure()
imagesc(inpar.xvalsFP,inpar.yvalsFP,PSFwBB);
axis image; 
axis([-3 3 -3 3]);
title('Broadband PSF w/ Vortex');
colorbar;
colormap(parula(256));
drawnow;