%{
Code for simulating the KPIC-VFN on-sky performance given known AO and other 
performance parameters. This script combines the results of the WF and PSF/EMap
scripts to generate the final SPIE figure.

This is based off the KPICVFN_OnSky_Simulator code.
    - This is the coupling analysis and figure generation part of that script.

Goals:___
- Simulate on-sky performance given synthetic or real WF residual data
- Account for atmospheric dispersion and ADC performance
- Generate figure

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
        ** [Updated]: Centering is based on star location since planet location
        is known only w.r.t star.
    * The eta_p reported in the bottom left graph is the average accross the
        full band based on the location that produced the best **average** (ie.
        broadband) planet coupling in the first frame.
        ** [Updated]: eta_p is calculated at user-provided location with
        user-provided radial bounds as before.
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


%-- Define Pupil phase (and ttres) files
% Folder with files
wfrfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\res12_apRad256_AllOn_118Samps_JuneRun\';
% filenames
wfrfnm = 'wfr';
ttrfnm = 'ttres';
% filename counter format
wfrFMT = '%06d';
% key params used in WFGenerator
wfStepSz = 500;
isrealWF = true;

%-- Define PSF and Coupling map files
% Folder with files
psfEMapfld = [wfrfld 'Charge2\'];
%filenames
psffnm = 'psfBB';    % psf files
eMapfnm = 'etaMaps';   % EMap files
% filename counter format
psfFMT = wfrFMT;
eMapFMT = wfrFMT;

%-- Define wfSamps based on number of fits files
wfSamps = size(dir([psfEMapfld eMapfnm '*.fits']),1);

% Flag to display SPIE figure or run in the background (true = display, false = hide)
isSPIEDisp = true;

%-- Saving parameters
% Flag to save gif
isSaveGif = true;
% Flag to save frames
isSaveFrames = true;
% Save folder
%svfld = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\KPIC_OnSkySims\IndividualScriptTest\';
svfld = [psfEMapfld 'Frames\'];
if isSaveGif
    % Create foler if it doesn't already exist
    mkdir(svfld)
end
% Save name for gif
svnmGif = 'SPIEFig.gif';
% Delay time for gif (time between frames in [s])
gifDelay = 0.3;
% Crop region for fits saving in [lambdaOverD]
fitcrop = 5;    

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

%% Define Figure properties
% Use subtightplot for control over plot spacing
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.06], [0.1 0.05], [0.05 0.05]);

%-- Define font and linewidth sizes
fontszAx = 14;
fontblAx = 'normal';
fontszTi = 14;
markerSz = 40;

%-- Axis limits
% PSF plot
axlimPSF = [-3 3 -3 3];     % x-y axis limits
% Eta Map
axlim2EMap = axlimPSF;      % x-y axis limits
% Eta_p plots
axYlimEtaP = [0 12];        % y axis limit
% Eta_s plots
isaxYlimEtaS = true;       % Flag to mark whether to hardcode y-axis
axYlimEtaS = [0 5];      % y axis limit (in units of eta_sYSCL)
eta_sYSCL  = 1e2;           % scaling for y axis
% realWF Zernike plot
zrnSCL = 1e2;               % scaling factor for y axis
axYlimZrn = [-10 10];          % y axis limits (in units of 1e3)

%-- Coloraxes
% WFR map
if isrealWF
    % Real WF thus get actual value directly from full matrix
        % slightly inaccurate since different mask and interp but close enough
    %caxlimWFR = [min(wfres(:)*lambda0*1e9) max(wfres(:)*lambda0*1e9)];
    caxlimWFR = [-400 400];
else
    % Approximate PV based on the RMS values from the coeffs
    caxlimWFR = [3*min(coeffs(:)*lambda0*1e9) 3*max(coeffs(:)*lambda0*1e9)];
end
% PSF
    % Determined empirically from several samples
    % Feel free to modify as needed
caxlimPSF = [0 0.22];  
% Eta map
caxlimEMap = axYlimEtaP/100;

%-- Central radius of planet coupling average region
eta_pR = 1.45;       %[lambda/D]

%-- Width of planet coupling average region
prg = 1;        % Average over region +/-prg pixels wide

%% Find optimal null in EMap files
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

%% Loop through samples
%-- Preallocate the eta_ss and eta_ps matrix for improved plotting
eta_ss = nan(wfSamps, 1);
eta_ps = nan(wfSamps, 1);

% Determine central band 
central_band_index = ceil(numWavelengths/2);

% Crop region for fits files
psfCrp = [N/2-round(lambdaOverD*5):N/2+round(lambdaOverD*5)];
etaMapCrp = [N/2-round(lambdaOverD*5):N/2+round(lambdaOverD*5)];

% Axes for plotting
xvals = -floor(dims(1)/2):floor(dims(1)/2)-1;
yvals = -floor(dims(2)/2):floor(dims(2)/2)-1;

% Size of cropped file (assuming dim1 - dim2
N2 = dims(1);

%-- Read ttres file
ttres = fitsread([wfrfld ttrfnm '.fits']);

for i = 1:wfSamps
fprintf('wfSamp %06d of %06d\n', i, wfSamps);
    
%%  Load data files

%-- Read WF File
wfr = fitsread([wfrfld sprintf([wfrfnm wfrFMT '.fits'],i)]);

%-- Read PSF file
iPSFv_BB = fitsread([psfEMapfld sprintf([psffnm psfFMT '.fits'],i)]);

%-- All eta_maps already read outside this loop; just extract frame
eta_maps = eMaps(:,:,:,i);

%% Do necessary anlyses
%-- Extract planet coupling (w.r.t star location)
eta_pL = nan(1,numWavelengths);
for ii = 1:numWavelengths
    %-- Azimuthally average at each wavelength, centering on star at each lam
    % Find star location (in pixels)
    cent = ttres(i, :, ii);     %[waves at lambda0]
    cent = cent*lambdaOverD;    %[pixels]
    % Average
    [tmpAvg, qvec] = VFN_An_radAverage(eta_maps(:,:,ii), round([N2/2 N2/2]+cent));
    % Find point closest to desired separation
    [~,ind2pull] = min(abs(qvec - eta_pR*lambdaOverD));
    eta_pL(ii) = mean(tmpAvg(ind2pull-prg:ind2pull+prg))*100;
end
eta_ps(i) = mean(eta_pL);

%-- Extract star coupling
eta_sL = squeeze(eta_maps(eta_sInd(1), eta_sInd(2), :))*eta_sYSCL;
eta_ss(i) = mean(eta_sL);

%% Make figure
if i == 1
    %-- First iteration: make the figure
    if isSPIEDisp
        gifFigVis = 'on';
    else
        gifFigVis = 'off';
    end

    gifFig = figure(7);
    set(gifFig, 'Units','normalized', 'OuterPosition', [0.05 0.05 0.9 0.9], 'Color', 'w', 'Visible', gifFigVis);
    
    %-- Plot instantaneous WFE/R (at central lambda)
    % Plot
    axWFR = subplot(2,3,1); 
    % Get temporary axes for wfr plot assuming equal dimensionality
    tmpax = size(wfr,1);
    tmpax = -tmpax/2:tmpax/2-1;
    datWFR = imagesc(tmpax/pupCircDiam,tmpax/pupCircDiam,wfr);
    axis image; 
    axis([-0.5 0.5 -0.5 0.5]);
    set(axWFR, 'TickDir', 'out');
    title(['Wavefront Residuals at ',num2str(lambdas(central_band_index)*1e9),'nm'], 'FontSize', fontszTi);
    xlabel('x/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    ylabel('y/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    wfrCBAR = colorbar; 
    wfrCBAR.Label.String = 'Residuals [nm]';
    wfrCBAR.Label.FontWeight = fontblAx;
    wfrCBAR.Label.FontSize = fontszAx;
    colormap(axWFR, parula(256));
    caxis(caxlimWFR);
    
    %-- Plot instantaneous Donut PSF (polychromatic)
    axPSF = subplot(2,3,2);
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
    colormap(axPSF, gray(256));
    caxis(caxlimPSF)
    
    %-- Plot eta_s vs. time
    axNull = subplot(2,3,5);
    if isrealWF
        % Rescale the x-axis to time (assuming 1/kHz sampling)
        xaxscl = wfStepSz/1000;
    else
        % Do not rescale
        xaxscl = 1;
    end
    datNull = scatter((1:wfSamps)*xaxscl, eta_ss, markerSz, 'b', 'MarkerFaceColor', 'flat');
    if ~isaxYlimEtaS
        axYlimEtaS(2) = 1.2*max(eta_sL(:));
    end
    ylim(axYlimEtaS)
    xlim([0 wfSamps*xaxscl])
    box on
    if eta_sYSCL == 1e2
        ylabel('\eta_s [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    else
        ylabel(['\eta_s [\times 10^{-' num2str(log10(eta_sYSCL)) '}]'], 'FontWeight', fontblAx, 'FontSize', fontszAx)
    end
    if isrealWF
        % Set x-axis label to seconds
        xlab = 'Time [s]';
    else
        xlab = 'WF Sample';
    end
    xlabel(xlab, 'FontWeight', fontblAx, 'FontSize', fontszAx)
    title('Broadband Star Coupling', 'FontSize', fontszTi)
    % Display running average (position empirically determined)
    eta_savg = annotation('textbox', [0.375 0.41 0.1 0.04], 'String', ['Average: ' num2str(mean(eta_ss,'omitnan'),'%0.2f') '%'], 'FontSize', fontszAx, 'FontWeight', fontblAx);

    %-- Plot eta_p vs time
    axPeak = subplot(2,3,4);
    datPeak = scatter((1:wfSamps)*xaxscl, eta_ps, markerSz, 'r', 'MarkerFaceColor', 'flat');
    ylim(axYlimEtaP)
    xlim([0 wfSamps*xaxscl])
    box on
    ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    xlabel(xlab, 'FontWeight', fontblAx, 'FontSize', fontszAx)
    title([sprintf('Broadband Planet Coupling\nfrom %3.1f', qvec(ind2pull-prg)/lambdaOverD) '\lambda/D to ' num2str(qvec(ind2pull+prg)/lambdaOverD,'%3.1f') '\lambda/D'], 'FontSize', fontszTi)
    % Display running average (position empirically determined) [0.055 0.11 0.105 0.04] <-- for 2-digit couplings
    eta_pavg = annotation('textbox', [0.055 0.11 0.10 0.04], 'String', ['Average: ' num2str(mean(eta_ps,'omitnan'),'%0.2f') '%'],'FontSize', fontszAx, 'FontWeight', fontblAx);
    
    %-- Plot instantaneous EMap (at lambda0)
    axEMap = subplot(2,3,3);
    datEMap = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,central_band_index));
    % Show circle where coupling was averaged. 
        % subtract -2 and xvals(end)... to center on image axes.
    crc1 = viscircles(([N2/2-1, N2/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (ind2pull-prg)/lambdaOverD, 'LineWidth', 1.0);
    set(crc1.Children(2), 'Visible', 'off');   % Hide outline of circle
    crc2 = viscircles(([N2/2-1, N2/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (ind2pull+prg)/lambdaOverD, 'LineWidth', 1.0);
    set(crc2.Children(2), 'Visible', 'off');   % Hide outline of circle
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
    colormap(axEMap, gray(256));
    caxis(caxlimEMap)
    
    %-- Plot eta_p and eta_s vs. wavelength
    axCoup = subplot(2,3,6);
    title('Coupling Fractions Across Band', 'FontSize', fontszTi)
    yyaxis left
    datCoupSLin = plot(lambdas*1e9, eta_sL, 'b-o', 'MarkerFaceColor', 'b');
    ylim(axYlimEtaS)
    if eta_sYSCL == 1e2
        ylabel('\eta_s [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    else
        ylabel(['\eta_s [\times 10^{-' num2str(log10(eta_sYSCL)) '}]'], 'FontWeight', fontblAx, 'FontSize', fontszAx)
    end
    xlabel('Wavelength [nm]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
    yyaxis right
    datCoupPLin = plot(lambdas*1e9, eta_pL, 'r-o', 'MarkerFaceColor', 'r');
    ylim(axYlimEtaP)
    ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
else
    %% Update figure
    %-- Other iterations:
    % Update WFR
    set(datWFR, 'CData', wfr);
    % Update Donut PSF
    set(datPSF, 'CData', iPSFv_BB);
    % Update eMap
    set(datEMap, 'CData', eta_maps(:,:,central_band_index));
    % Update circles on eMap
    crc1.delete; crc2.delete;
    crc1 = viscircles(axEMap,(([N2/2-1, N2/2-1]+cent)/lambdaOverD-xvals(end)/lambdaOverD),  (ind2pull-prg)/lambdaOverD, 'LineWidth', 1.0);
    set(crc1.Children(2), 'Visible', 'off');   % Hide outline of circle
    crc2 = viscircles(axEMap,(([N2/2-1, N2/2-1]+cent)/lambdaOverD-xvals(end)/lambdaOverD),  (ind2pull+prg)/lambdaOverD, 'LineWidth', 1.0);
    set(crc2.Children(2), 'Visible', 'off');   % Hide outline of circle
    % Update eta_s vs. time
    set(datNull, 'YData', eta_ss);
    eta_savg.String = ['Average: ' num2str(mean(eta_ss,'omitnan'),'%0.2f') '%'];
    % Update eta_p vs. time
    set(datPeak, 'YData', eta_ps);
    eta_pavg.String = ['Average: ' num2str(mean(eta_ps,'omitnan'),'%0.2f') '%'];
    % Update etas vs. wavelength
    set(datCoupSLin, 'YData', eta_sL);
    set(datCoupPLin, 'YData', eta_pL);
end

drawnow;

if isSaveGif
    % Save initial frame
    frame = getframe(gifFig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1
        imwrite(imind,cm,[svfld svnmGif],'gif', 'Loopcount',inf, 'DelayTime', gifDelay);
    else
        imwrite(imind,cm,[svfld svnmGif],'gif','WriteMode','append','DelayTime', gifDelay); 
    end
    if isSaveFrames
        saveas(gifFig, [svfld sprintf(['Frame' wfrFMT '.png'],i)]);
    end
end

end

if isSaveGif
    saveas(gifFig, [svfld 'SPIEFig_FinalFrame.png'])
end