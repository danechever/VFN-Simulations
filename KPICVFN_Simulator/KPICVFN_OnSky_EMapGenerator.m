%{
Code for simulating the KPIC-VFN on-sky performance given known AO and other 
performance parameters. This script produces the PSF and Coupling maps.

This is based off the KPICVFN_OnSky_Simulator code.
    - This is the analysis section of that code.

Goals:___
- Simulate on-sky performance given synthetic or real WF residual data
- Account for atmospheric dispersion and ADC performance

Notes:___
* Wavefront residuals are plotted in [nm] at central lambda
* Eta_s is calculated using min in central region:
    * Multi-wavelength eta matrix gets cropped to central subsection of frame
    * All wavelengths are averaged to get broadband null 
    * Min of this wavelength-averaged sub-matrix is treated as null
    --> Thus, the calculated null puts the fiber at the best location for a
        **broadband** null assuming a flat spectrum. 
    --> This has implications to be considered when the null splits.
    * The wavelength eta_s's (for the final plot) are calculated with the fiber
        at the location that produced the best broadband null, NOT the ideal
        position for each individual wavelength.
* Y-axis of eta_s plots are rescaled by a user-provided value to avoid having
    the extra exponent at the top of the figure.
    * [Updated]: User can also select to have axes dynamically scaled by code to
    go from 0 to 1.2*max(eta_sL).
* Eta_p is calculated via azimuthal averaging
    * VFN_An_radAverage (from VFN-Lab repo) is used for this averaging
    * The centering of the average is based on the **broadband** null location
        calculated as described above. 
        ** [Updated]: Centering is based on center of frame since this is where
        the star is actually located and hence the planet distance is w.r.t here
    * The eta_p reported in the bottom left graph is the average accross the
        full band based on the location that produced the best **average** (ie.
        broadband) planet coupling in the first frame.
    * The wavelength-dependent couplings (for the final plot) are calculated
        w.r.t the real star location as described in last bullet but at each 
        wavelength.
    * [Updated]: Averaging is now not just at a single separation, it is over a
        user-defined range of radii. This allows for uncertainty in the planet
        location.
* A circle is superimposed on the eta map to show where the ideal planet
    location is based on the Eta_p calculation described above. The average at
    this circle's location is the value reported as eta_p everywhere else.
    ** [Updated]: Two circles are now used to show the inner and outer range of
    the averaged region.
* Keck pupil is used instead of simple circle.
* Synthetic WF residuals are modulated accross samples using sine function to 
    have continuously variable WF.
    ** [Updated] This is now optional.
* Magnitude of each WF residual across samples is plotted at the end.
    ** [Updated] This is only done when a synthetic WF is used
    ** [Updated] This is done for both synthetic and real wavefronts and is done
        at the beginning (for synthetic) or live during the loop (for real)
* "subtightplot" library is used to provide control of subplot location.
    * This is NOT a native matlab library; user must download from file exchange
* Full figure is only plotted once and then only the data is updated on each
    iteration of wavefront samples.
* Gif is optionally output at the end.
* For real WF data:
    - User must provide fits files with [xdim, ydim, sampleDim] format
    - User also provides a course mask for this data
        - That mask is only used for identifying which pixels/actuators are
        actually within the real beam.
    - The WF data is rotated and rescaled to match the sampling of our
        simulation. This is done using interpolation. I've tried to do these
        steps in such a way that the data is minimally affected.
    - Once the WF is oriented and scaled correctly, our simulated keck pupil
        (which includes central obscuration, spiders, fine outer pattern, etc.)
        is applied to the WF data. 
    ** [Updated] Fudge-factors have been added for user to tune the centering,
        orientation, and size of the simulated pupil within the WF data to
        ensure it is sampling the same region as the real pupil mask.
    ** [Updated] the raw WF (described above) can be used directly or a zernike
        decomp can be done and this reconstructed wavefront can be applied
        instead.
* Plotting is now done much more efficiently:
    - As before, main figure is only produced once and then updated.
    - However, instead of saving all frames of all the data, only the current
        frame's is saved and this gets overwritten on each WF sample iteration.
    - This vastly decreases RAM use and increases calculation speed.
* Code now handles tip/tilt residuals as well as ADC residulas
    - User can input a tip/tilt residual fits file to process real data. This
        file must be in [arcsec] and have dimensionality [wfSamps, ti[/tilt]
    - Synthetic residuals can also be provided in [mas] directly. These can be
        modulated in the code or kept constant.
    - ADC residuals are provided at each sample and wavelength.
%}

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
charge = 1*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

%-- Define Pupil phase files
% Folder with files
phzfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\res12_apRad256_AllOn_118Samps_JuneRun\';
% filename
phzfnm = 'pupphz';
% filename counter format
phzFMT = '%06d';

%-- Define wfSamps based on number of fits files
wfSamps = size(dir([phzfld phzfnm '*.fits']),1);

%--Offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = true;

%-- Saving parameters
% Flag to save gif
isSaveGif = true;
% Flag to save fits
isSaveFit = true;
% Save folder
%svfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\IndividualScriptTest\';
svfld = [phzfld 'Charge1\'];
if isSaveFit
    % Create foler if it doesn't already exist
    mkdir(svfld)
end
% Save name for gif
svnmGif = 'PSF_EMaps.gif';
% Delay time for gif (time between frames in [s])
gifDelay = 0.3;
% Crop region for fits saving in [lambdaOverD]
fitcrop = 5;    

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/(pupCircDiam/2)/2; % lam/D in units of samples in the image plane

% figure(1)
% imagesc(xvals/(pupCircDiam/2),yvals/(pupCircDiam/2),PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Make vortex mask (once before loop since does not change)
EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

%% Generate fibermode at each lambda (once before loop since does not change)
    % Generating fib modes once improves runtime by 25% for N=2048 and wfSamps=10
    
% Parameters for Thorlabs SM600
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

%% Loop through WF samples
% Use subtightplot for control over plot spacing
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.06], [0.1 0.05], [0.05 0.05]);

%-- Define font and linewidth sizes
fontszAx = 14;
fontblAx = 'normal';
fontszTi = 14;

%-- Axis limits
% PSF plot
axlimPSF = [-3 3 -3 3];     % x-y axis limits
% Eta Map
axlim2EMap = axlimPSF;      % x-y axis limits
% Eta_p plots
axYlimEtaP = [0 20];        % y axis limit

%-- Coloraxes
% PSF
    % Determined empirically from several samples
    % Feel free to modify as needed
caxlimPSF = [0 0.33];  
% Eta map
caxlimEMap = axYlimEtaP/100;

% Determine central band 
central_band_index = ceil(numWavelengths/2);

% Crop region for fits files
psfCrp = [N/2-round(lambdaOverD*5):N/2+round(lambdaOverD*5)];
etaMapCrp = [N/2-round(lambdaOverD*5):N/2+round(lambdaOverD*5)];

for i = 1:wfSamps
fprintf('wfSamp %06d of %06d\n', i, wfSamps);
    
%%  Define pupil field

%-- Read WF File
pupphz = fitsread([phzfld sprintf([phzfnm phzFMT '.fits'],i)]);

%-- Pad to match PUPIL size
% Pad to match Simulated pupil diameter
pd = size(pupphz(:,:,1)) - [N N];
if pd(1) < 0
    % Need to add rows 
    pupphz = padarray(pupphz, [-round(pd(1)/2) 0 0]);
end
if pd(2) < 0
    % Need to add columns
    pupphz = padarray(pupphz, [0 -round(pd(2)/2) 0]);
end
% Final crop to ensure proper size: padding adds 1 extra row/col
pupphz = pupphz(1:N, 1:N,:);

for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*pupphz(:,:,ch)).*PUPIL;
end

%% Calculate PSF

iPSFv_BB = getPSF(Epup.*EPM,lambda0,lambdas,normI,coords);

%% Generate coupling maps
% Use pre-made fibermodes to improve runtime
eta_maps = generateCouplingMap_polychromatic(Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords, fibmodes);

%% Plot if desired
if isPlotSimp
    %-- Display PSF (Broadband) and EMap (lambda0) (top row in SPIE figure)
    if i == 1
        %-- First iteration: make the figure
        fig = figure(7);
        set(fig, 'units', 'inches', 'Position', [0 0 14 6], 'color', 'w');
        %-- Plot instantaneous Donut PSF (polychromatic)
        axPSF = subplot(1,2,1);
        datPSF = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
        axis image; 
        set(axPSF, 'TickDir', 'out');
        axis(axlimPSF);
        title('Broadband PSF w/ Vortex', 'FontSize', fontszTi);
        xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        psfCBAR = colorbar; 
        psfCBAR.Label.String = 'Normalized Irradiance';
        psfCBAR.Label.FontWeight = fontblAx;
        psfCBAR.Label.FontSize = fontszAx;
        %caxis([-3 0])
        %set(axPSF, 'Colormap', gray(256))
        colormap(axPSF, gray(256));
        caxis(caxlimPSF)
        
        %-- Plot instantaneous EMap (at lambda0)
        axEMap = subplot(1,2,2);
        datEMap = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,central_band_index));
        % Can't show circle where eta_p will be calculated since line-profiles
        % are not generated in this script, thus no peak is calculated
        axis image;
        axis(axlim2EMap);
        set(axEMap, 'TickDir', 'out');
        title(['Coupling at ',num2str(lambda0*1e9),'nm'], 'FontSize', fontszTi);
        xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
        EMapCBAR = colorbar; 
        EMapCBAR.Label.String = '\eta';
        EMapCBAR.Label.FontWeight = fontblAx;
        EMapCBAR.Label.FontSize = fontszAx;
        %set(axEMap, 'Colormap', gray(256))
        colormap(axEMap, gray(256));
        caxis(caxlimEMap)
    else
        %-- Other iterations:
        % Update Donut PSF
        set(datPSF, 'CData', iPSFv_BB);
        set(datEMap, 'CData', eta_maps(:,:,central_band_index));
    end

    drawnow;
end
if isSaveFit
    %-- Save PSF
    fitswrite(iPSFv_BB(psfCrp,psfCrp), [svfld sprintf('psfBB%06d.fits',i)]);
    %-- Save EMaps
    fitswrite(eta_maps(etaMapCrp,etaMapCrp,:), [svfld sprintf('etaMaps%06d.fits',i)]);
end

if isSaveGif && isPlotSimp
    % Save initial frame
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1
        imwrite(imind,cm,[svfld svnmGif],'gif', 'Loopcount',inf, 'DelayTime', gifDelay);
    else
        imwrite(imind,cm,[svfld svnmGif],'gif','WriteMode','append','DelayTime', gifDelay); 
    end
end

end

if isSaveGif && isPlotSimp
    saveas(fig, [svfld 'PSF_Emaps_FinalFrame.png'])
end