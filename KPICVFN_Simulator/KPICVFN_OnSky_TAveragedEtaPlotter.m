clear; close all; 
addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))

%% Input parameters 

% Define smapling info
N = 2^12; % Size of computational grid (NxN samples) 
apRad = 128*2; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2200e-9; %central wavelength
fracBW = 0.20; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

%--Offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

%-- Define PSF and Coupling map files
% Folder with files
wfrfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\res12_apRad256_AllOn_118Samps_JuneRun\';
psfEMapfld = [wfrfld 'Charge2\'];
%filenames
psffnm = 'psfBB';    % psf files
eMapfnm = 'etaMaps';   % EMap files
% filename counter format
psfFMT = '%06d';
eMapFMT = psfFMT;

%-- Saving parameters
% Flag to save plot
isSaveFig = true;
% Save folder
svfld = 'C:\Users\danie\Desktop\';
% Save name for plot
svnmFig = 'Charge2_TAveragedEMapTEST.png';
 

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function
disp('Making PUPIL')
%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/(pupCircDiam/2)/2; % lam/D in units of samples in the image plane

%% Make vortex mask 
disp('Generating Vortex')
EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

%% Generate fibermode at each lambda
disp('Making Fiber Modes')    

% Parameters for Thorlabs SM2000
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4436;% core index 
fiber_props.n_clad = 1.4381;% cladding index 
Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.4); % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

% Iterate through wavelengths generating modes
fibmodes = nan(N, N, numWavelengths);
for ch = 1:numWavelengths
    fibmodes(:,:,ch) = generateSMFmode(fiber_props, lambdas(ch), Fnum, lambdaOverD, coords);
end

%% Generate ideal coupling maps
%-- Def pupil phase
for ch = 1:numWavelengths
    Epup(:,:,ch) = PUPIL;
end

disp('Calculating Ideal Coup Maps')
%-- Get average coupling maps
eMapIdeal = generateCouplingMap_polychromatic(Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords, fibmodes);
eMapIdeal = mean(eMapIdeal, 3);

%% Get simulated average (coupling maps should have already been generated)
disp('Loading PyWFS Coup Maps')

%-- Define wfSamps based on number of fits files
wfSamps = size(dir([psfEMapfld eMapfnm '*.fits']),1);

%-- Pre-allocate matrix for EMap files
% Open single file to get dimensions
dims = fitsinfo([psfEMapfld eMapfnm sprintf([eMapFMT '.fits'],1)]);
dims = cell2mat(dims.PrimaryData.Keywords(4:6,2))';
eMaps = nan([dims wfSamps]);

%-- Load all EMap files
for i = 1:wfSamps
    eMaps(:,:,:,i) = fitsread([psfEMapfld eMapfnm sprintf([eMapFMT '.fits'],i)]);
end

%-- Average along spectral and temporal dimensions to get optimal fiber position
eMapmean = mean(mean(eMaps, 3), 4);

% Get eta_s as min in center of all eta_maps averaged together
    % Average eta_maps since fiber can only be in one location so need to make
    % sure we take the eta_s at this same locaiton always.
crp = round(lambdaOverD);   % Use 1lam/D as 1/2 of crop window; 
                   %assumes null stays within 1lam/D on average
[s1, s2] = size(eMapmean);
eta_cent = eMapmean(ceil(s1/2-s1/crp):ceil(s1/2+s1/crp),ceil(s2/2-s2/crp):ceil(s2/2+s2/crp));
eta_s = min(eta_cent(:));    % Get min of average in central region as opt. null


%-- Get indices for optimal null location
[eta_sInd(1), eta_sInd(2)] = find(eta_s==eta_cent);
eta_sInd = eta_sInd + [floor(s1*(1/2-1/crp)) floor(s2*(1/2-1/crp))];   % Account for cropping 

%% Plotting Parameters
disp('Plotting')
% Use subtightplot for control over plot spacing
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.06], [0.1 0.05], [0.05 0.05]);

%-- Define font and linewidth sizes
fontszAx = 14;
fontblAx = 'normal';

%-- FOV for plots
fovReg = [-3 3 -3 3];
fovLog = [-5 5 -5 5];


%% Plot
%-- Create figure
fig = figure('Color', 'white', 'Units', 'Inches', 'Position', [0 0 10 10]);

%-- SIMULATED --
% Set the axes size
xvals = -floor(dims(1)/2):floor(dims(1)/2)-1;
yvals = -floor(dims(2)/2):floor(dims(2)/2)-1;

%- Regular
sp1 = subplot(2, 2, 1);
ax1 = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eMapmean);
axis image; colormap(sp1, gray(256));
title('Time-Averaged \eta Map')
axis(fovReg);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
cbr1 = colorbar;
set(cbr1.Label, 'String', '\eta', 'FontWeight', fontblAx, 'FontSize', fontszAx);

%- Log10
sp2 = subplot(2, 2, 2);
ax2 = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eMapmean));
axis image; colormap(sp2, gray(256));
title('Time-Averaged \eta Map (log10)')
axis(fovLog);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
cbr2 = colorbar;
set(cbr2.Label, 'String', 'log10(\eta)', 'FontWeight', fontblAx, 'FontSize', fontszAx);

%-- IDEAL --
% Set the axes size
xvals = coords.xvals;
yvals = coords.yvals;

%- Regular
sp3 = subplot(2, 2, 3);
ax3 = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eMapIdeal);
axis image; colormap(sp3, gray(256));
title('Time-Averaged Ideal \eta Map')
axis(fovReg);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
cbr3 = colorbar;
caxis(sp3, sp1.CLim)
set(cbr3.Label, 'String', '\eta', 'FontWeight', fontblAx, 'FontSize', fontszAx);

%- Log10
sp4 = subplot(2, 2, 4);
ax4 = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eMapIdeal));
axis image; colormap(sp4, gray(256));
title('Time-Averaged Ideal \eta Map (log10)')
axis(fovLog);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx);
cbr4 = colorbar;
caxis(sp4, sp2.CLim)
set(cbr4.Label, 'String', 'log10(\eta)', 'FontWeight', fontblAx, 'FontSize', fontszAx);

%% Save as needed
if isSaveFig
    export_fig([svfld svnmFig], '-r300', '-painters');
end