%{
Code for simulating the KPIC-VFN on-sky performance given known AO and other 
performance parameters. This code also produces the SPIE 2019 figure.

This is based off the PolychromaticSimTests code.

Goals:___
- Simulate on-sky performance given synthetic or real WF residual data
- Account for atmospheric dispersion and ADC performance
- Efficiently produce final figure

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
* Magnitude of each WF residual across samples is plotted at the end.
    ** [Updated] This is only done when a synthetic WF is used
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
* Plotting is now done much more efficiently:
    - As before, main figure is only produced once and then updated.
    - However, instead of saving all frames of all the data, only the current
        frame's is saved and this gets overwritten on each WF sample iteration.
    - This vastly decreases RAM use and increases calculation speed.
%}

clear; %close all; 
addpath('VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))

%% Input parameters 

% Define smapling info
N = 2^11; % Size of computational grid (NxN samples) 
apRad = 128; % Aperture radius in samples 

% Define wavelength info
lambda0 = 2200e-9; %central wavelength
fracBW = 0.20; %\Delta\lambda/\lambda
numWavelengths = 3;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

%--- Define wavefront error at the central wavelength
isrealWF = false;      % Flag to mark if real or synthetic WFs should be used
if isrealWF
    % Real WF data parameters
    wfPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
    wfRSN = 'residualWF.fits';          % Filename for residuals
    wfMKN = 'DMPupil.fits';             % Filename for pupil mask on residuals
    % Define sample start and end values
    wfstrt = 100;
    wfSamps = 20;
    % Scaling factor for data (nm/V)
    wfSCL = 600;     
    %%%% NOTE:::: DEAL WITH MANUAL CORRECTIONS IN WF SECTION BELOW!!!
else
    % Synthetic wavefront parameters
    nolls = [2, 5, 6];
    coeffs = 1*[0.01, 0.025, 0.01];    % Waves RMS
    wfSamps = 10;                      % Number of WF samples to synthesize
    wfModPhzEnd = 2*pi;                % Phase end (without delay)
end

% Give offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

% Flag to plot coupling maps in log scale
islogcoup = true;

% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = false;

%-- Saving parameters
% Flag to save data
isSaveGif = false;
% Save folder
svfld = 'C:\Users\danie\Desktop\';
svnm  = 'RealWF_noTT_Samp100-120.gif';

%% Generate the coordinate system

coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%% Create array with pupil function

%PUPIL = makeCircularPupil( apRad, N );
[PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );
[normI, totalPower0] = getNormalization(PUPIL);% Normalization factors
lambdaOverD = N/apRad/2; % lam/D in units of samples in the image plane

% figure(1)
% imagesc(xvals/apRad,yvals/apRad,PUPIL);
% axis image; 
% axis([-1 1 -1 1]);
% title('Pupil');
% colorbar; 
% colormap(parula(256));
% drawnow;

%% Deal with WF extraction (real data) or synthesis (synthetic data)
if isrealWF
    %*** Real data
    
    %-- Read data 
    wfres = fitsread([wfPTH, wfRSN]);   % Full residuals file
    wfmsk = fitsread([wfPTH, wfMKN]);   % Mask file
    % Extract specific samples
    wfres = wfres(:,:,wfstrt:wfstrt+wfSamps-1);
    % Rescale the WFR to waves at lambda0
    wfres = wfres*wfSCL*1e-9/lambda0;    %x1e-9 to go to [m] since lambda0 is in [m]
    % Copy raw mask for plotting later
    wfmskRAW = wfmsk;
    % Extract in-pupil pixels (raw file is 0=out, 1=in pupil, 2=in DM)
    wfmsk(wfmsk~=1) = 0;
    
    %-- Mask corrections (apply manually on case-by-case basis)
    % Fill in central pixels
    wfmsk(10:12, 9:10) = 1;
    % Correct shifted pixels (1 extra at top and one missing at bottom)
    wfmsk(2,13) = 0; wfmsk(19,13) = 1;
    
    %-- Pre-process WF before loop
    % Mask WF using charlotte's pupil
    wfmsk = repmat(wfmsk, [1 1 wfSamps]);
    wfres = wfres.*wfmsk;
    % Pad sides of WF so that up-sampling produces even-sized grid
    wfres = padarray(wfres, [1 1 0]);
    
else
    %*** Synthetic data
    % Create matrix-version of coeffs to simulate changing WFE
    wfMod = nan(wfSamps, length(nolls));
    for i = 1:length(nolls)
        % Create random phase delay
        phzShift = rand*2*pi;
        wfMod(:, i) = sin(linspace(0+phzShift, wfModPhzEnd+phzShift, wfSamps))';    % smoothly modulate WFE amplitude
    end
    %wfMod = repmat(wfMod, 1, length(coeffs));
    coeffs = repmat(coeffs, wfSamps, 1);
    coeffs = coeffs.*wfMod;             % Apply modulation
end

%% Display WF Mask (for real data) or Zernikes (for synthetic data)
if isrealWF
    if isPlotSimp
        % Show raw mask
        figure(100);
        imagesc(wfmskRAW);
        axis image;
        % Show corrected mask
        figure(101);
        imagesc(wfmsk(:,:,1))
        axis image;
    end
else
    % Show synthetic coefficients
    zrnFig = figure(100);
    mrks = {'-o', '-s', '-d', '-^', '-*', '-x', '-p'};
    legs = cell(length(nolls),1);
    % plot(1:wfSamps, coeffs, mrks)
    for i = 1:length(nolls)
        plot(coeffs(:,i), char(mrks(i)), 'MarkerSize', 4, 'LineWidth', 2);
        if i == 1
            hold on
        end
        legs(i) = {num2str(nolls(i))};
    end
    hold off

    legH = legend(legs);
    title(legH, 'Noll Index')
end    

%% Make vortex mask (once before loop since does not change)
EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

%% Set SPIE Figure parameters
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
axYlimEtaS = [0 1];      % y axis limit (in units of 1e-3)
eta_sYSCL  = 1e3;           % scaling for y axis

%-- Coloraxes
% WFR map
if isrealWF
    % Real WF thus get actual value directly from full matrix
        % slightly inaccurate since different mask and interp but close enough
    caxlimWFR = [min(wfres*lambda0*1e9, [], 'all') max(wfres*lambda0*1e9, [], 'all')];
else
    % Approximate PV based on the RMS values from the coeffs
    caxlimWFR = [3*min(coeffs*lambda0*1e9,[],'all') 3*max(coeffs*lambda0*1e9, [], 'all')];
end
% PSF
    % Determined empirically from several samples
    % Feel free to modify as needed
caxlimPSF = [0 0.22];  
% Eta map
caxlimEMap = axYlimEtaP/100;

%-- Width of planet coupling average region
prg = 1;        % Average over region +/-prg pixels wide

%% Loop through WF samples
%-- Preallocate the eta_ss and eta_ps matrix for improved plotting
eta_ss = nan(wfSamps, 1);
eta_ps = nan(wfSamps, 1);

for i = 1:wfSamps
fprintf('wfSamp %03d of %03d\n', i, wfSamps);
    
%% Define pupil field
if isrealWF   
    % Extract frame
    wfTMP = wfres(:,:,i);
    % Up-sample to decrease error sensitivity in centering and rotation
        % Want enough pts that rot. interp. via nearest-neighbor looks good.
    wfTMP = interp2(wfTMP, 5, 'nearest');
    % Recenter and crop (assuming roughly symmetric)
    [rws, cls] = find(wfTMP);  % find min and max row/col indices
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1-1:rw2, cl1-1:cl2);
    % Rotate to fix orientation
    wfTMP = imrotate(wfTMP, 30);
    % Find inscribed circle size
    [rws, cls] = find(~wfTMP);
    crc = 2*min(sqrt((rws-size(wfTMP,1)/2-1).^2 + (cls-size(wfTMP,2)/2-1).^2));
    % OPTIONAL: show inscribed circle
    %figure; imagesc(wfTMP)
    %viscircles([size(wfTMP,2)/2-1, size(wfTMP,1)/2-1], crc/2);
    % Resize so that inscribed circle is >~= pupCircDiam
    fdg = 5;    % fudge factor on pupCircDiam to make sure bigger than crc
	x = 0:(pupCircDiam+fdg)/crc:(size(wfTMP,2)-1)*(pupCircDiam+fdg)/crc;
    y = 0:(pupCircDiam+fdg)/crc:(size(wfTMP,1)-1)*(pupCircDiam+fdg)/crc;
    xi = 0:(size(wfTMP,2)-1)*(pupCircDiam+fdg)/crc;
    yi = 0:(size(wfTMP,1)-1)*(pupCircDiam+fdg)/crc;
    [xt, yt] = meshgrid(x,y); %create vector arrays
    [xit, yit] = meshgrid(xi,yi); %create vector arrays
    wfTMP = interp2(xt,yt,wfTMP,xit,yit);
    % Pad to match Simulated pupil diameter
    pd = size(wfTMP) - [N N];
    if pd(1) < 0
        % Need to add rows 
        wfTMP = padarray(wfTMP, [-round(pd(1)/2) 0]);
    end
    if pd(2) < 0
        % Need to add columns
        wfTMP = padarray(wfTMP, [0 -round(pd(2)/2)]);
    end
    % Final crop to ensure proper size 
        % Accounts for initial oversize or padding which adds 1 extra row/col
    wfTMP = wfTMP(1:N, 1:N);
    % Apply pupil!  (Not actually, pupil will be applied at Epup calc
    % wfTMP = wfTMP.*PUPIL;
    % OPTIONAL: show final pupil with WF
    %figure; imagesc(wfTMP.*PUPIL)
    
    % Convert from waves to radians
    phz = wfTMP*2*pi;
else
    phz = generateZernike_fromList( nolls, coeffs(i,:), PUPIL, pupCircDiam/2, coords); 
end

if isPlotSimp
    figure(2);
end
for ch = 1:numWavelengths
    Epup(:,:,ch) = exp(1i*phz*lambda0/lambdas(ch)).*PUPIL;
    
    if isPlotSimp
        subplot(1,numWavelengths,ch);
        imagesc(xvals/apRad,yvals/apRad,angle(Epup(:,:,ch)));
        axis image; 
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(parula(256));
    end
drawnow;   
end

%% Get PSF without vortex mask

% Get broadband PSF
iPSF_BB = getPSF(Epup,lambda0,lambdas,normI,coords);

if isPlotSimp
    figure(3)
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSF_BB);
    axis image; 
    axis([-3 3 -3 3]);
    title('broadband PSF w/o vortex');
    colorbar;%caxis([-3 0])
    colormap(parula(256));
    drawnow;
end 

%% Plot vortex mask
central_band_index = ceil(numWavelengths/2);

if isPlotSimp
    figure(4)
    imagesc(xvals/apRad,yvals/apRad,angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
    axis image; 
    axis([-1 1 -1 1]);
    title('Pupil phase at \lambda_0');
    colormap(hsv(256));
    colorbar; 
    drawnow;
end

%% Get PSF with vortex mask

iPSFv_BB = getPSF(Epup.*EPM,lambda0,lambdas,normI,coords);

if isPlotSimp
    figure(5)
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
    axis image; 
    axis([-3 3 -3 3]);
    title('broadband PSF w/ vortex');
    colorbar;%caxis([-3 0])
    colormap(parula(256));
    drawnow;
end

%% Generate coupling maps for each wavelength

% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
Fnum = 5; % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

eta_maps = generateCouplingMap_polychromatic(Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords);

if isPlotSimp
    figure(6);
    for ch = 1:numWavelengths
        subplot(1,numWavelengths,ch);
        if islogcoup
            imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eta_maps(:,:,ch)));
        else
            imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,ch));
        end
        axis image;
        axis([-2 2 -2 2]);
        title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
        colorbar;
        if islogcoup; caxis([-4 0]); end
        colormap(gray(256));
    end
drawnow
end


%% Deal with SPIE figure

if i == 1
    %% Generate figure (first sample only)
gifFig = figure(7);
set(gifFig, 'Units','normalized', 'OuterPosition', [0.05 0.05 0.9 0.9], 'Color', 'w');

%-- Plot instantaneous WFE/R (at central lambda)
% Convert WFE from rads to nm
wfr = phz/2/pi*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
% Set out-of-pupil values to NaN for plotting
wfr(wfr==0) = nan;
% Plot
axWFR = subplot(2,3,1); 
datWFR = imagesc(xvals/pupCircDiam,yvals/pupCircDiam,wfr);
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
set(axWFR, 'Colormap', parula(256))
caxis(caxlimWFR)

%-- Plot instantaneous Donut PSF (polychromatic)
axPSF = subplot(2,3,2);
datPSF = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFv_BB);
axis image; 
set(axPSF, 'TickDir', 'out');
axis(axlimPSF);
title('Broadband PSF w/ vortex', 'FontSize', fontszTi);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
psfCBAR = colorbar; 
psfCBAR.Label.String = 'Normalized Irradiance';
psfCBAR.Label.FontWeight = fontblAx;
psfCBAR.Label.FontSize = fontszAx;
%caxis([-3 0])
set(axPSF, 'Colormap', gray(256))
caxis(caxlimPSF)

%-- Plot eta_s vs. time [eta_s first since it's needed for eta_p and eta map]
% Get eta_s as min in center of all eta_maps averaged together
    % Average eta_maps since fiber can only be in one location so need to make
    % sure we take the eta_s at this same locaiton always.
crp = 68;   %1/2 of crop window for min analysis
[s1, s2, ~] = size(eta_maps);
eta_cent = eta_maps(ceil(s1/2-s1/crp):ceil(s1/2+s1/crp),ceil(s2/2-s2/crp):ceil(s2/2+s2/crp),:);
eta_cent = mean(eta_cent,3);
eta_s = min(eta_cent,[],[1 2]);    % Get min of average in central region as broadband null
% Get indices for null location
[eta_sInd(1), eta_sInd(2)] = find(eta_s==eta_cent);
eta_sInd = eta_sInd + round(N/2-N/crp-1);   % Account for cropping from before
% Eta_s vs wavelength (at eta_sInd location)
eta_sL = squeeze(eta_maps(eta_sInd(1), eta_sInd(2), :))*eta_sYSCL;
% Get broadband null value (via average)
eta_ss(1) = mean(eta_sL);   
% Plot
axNull = subplot(2,3,5);
datNull = scatter(1:wfSamps, eta_ss, markerSz, 'b', 'MarkerFaceColor', 'flat');
if ~isaxYlimEtaS
    axYlimEtaS(2) = 1.2*max(eta_sL(:));
end
ylim(axYlimEtaS)
xlim([0 wfSamps])
box on
ylabel('\eta_s [\times 10^{-3}]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('WF Sample', 'FontWeight', fontblAx, 'FontSize', fontszAx)
title('Broadband Star Coupling', 'FontSize', fontszTi)

%-- Plot eta_p vs. time [eta_p before eta map since needed for circle plot]
% Iterate through wavelengths to get azimuthal average of coupling
%eta_sI = nan(numWavelengths, 2);
eta_pAvg = nan(numWavelengths, N/2-1);
for ii = 1:numWavelengths
    % Azimuthally average
    [tmpAvg, qvec] = VFN_An_radAverage(eta_maps(:,:,ii), [N/2+1, N/2+1]);
    eta_pAvg(ii,:) = tmpAvg;
end
% Average along wavelength dimension to find best average planet coupling
[~, eta_pInd] = max(mean(eta_pAvg, 1));
eta_ps(1) = mean(eta_pAvg(:,eta_pInd-prg:eta_pInd+prg), [1 2])*100;  % convert to [%]
% Plot
axPeak = subplot(2,3,4);
datPeak = scatter(1:wfSamps, eta_ps, markerSz, 'r', 'MarkerFaceColor', 'flat');
ylim(axYlimEtaP)
xlim([0 wfSamps])
box on
ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('WF Sample', 'FontWeight', fontblAx, 'FontSize', fontszAx)
title([sprintf('Broadband Planet Coupling\nfrom %3.1f', qvec(eta_pInd-prg)/lambdaOverD) '\lambda/D to ' num2str(qvec(eta_pInd+prg)/lambdaOverD,'%3.1f') '\lambda/D'], 'FontSize', fontszTi)

%-- Plot instantaneous eta map (at central lambda)
axEMap = subplot(2,3,3);
datEMap = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,central_band_index));
% Show circle where coupling was averaged. 
    % subtract -2 and xvals(end)... to center on image axes.
crc1 = viscircles(([N/2-1, N/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (eta_pInd-prg)/lambdaOverD, 'LineWidth', 1.0);
set(crc1.Children(2), 'Visible', 'off');   % Hide outline of circle
crc2 = viscircles(([N/2-1, N/2-1]/lambdaOverD-xvals(end)/lambdaOverD),  (eta_pInd+prg)/lambdaOverD, 'LineWidth', 1.0);
set(crc2.Children(2), 'Visible', 'off');   % Hide outline of circle
axis image;
axis(axlim2EMap);
set(axEMap, 'TickDir', 'out');
title(['Coupling at ',num2str(lambdas(central_band_index)*1e9),'nm'], 'FontSize', fontszTi);
xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
EMapCBAR = colorbar; 
EMapCBAR.Label.String = '\eta';
EMapCBAR.Label.FontWeight = fontblAx;
EMapCBAR.Label.FontSize = fontszAx;
set(axEMap, 'Colormap', gray(256))
caxis(caxlimEMap)

%-- Plot eta_p and eta_s vs. wavelength
% Eta_p vs wavelength (average w/ centering based on eta_sInd)
eta_pL = mean(eta_pAvg(:,eta_pInd-prg:eta_pInd+prg),2)*100;
% Plot
axCoup = subplot(2,3,6);
title('Coupling Fractions Across Band', 'FontSize', fontszTi)
yyaxis left
datCoupSLin = plot(lambdas*1e9, eta_sL, 'b-o', 'MarkerFaceColor', 'b');
%hold on
%datCoupSSca = scatter(lambdas*1e9, eta_sLs(1,:));
%hold off
ylim(axYlimEtaS)
ylabel(['\eta_s [\times 10^{-' num2str(log10(eta_sYSCL)) '}]'], 'FontWeight', fontblAx, 'FontSize', fontszAx)
xlabel('Wavelength [nm]', 'FontWeight', fontblAx, 'FontSize', fontszAx)
yyaxis right
datCoupPLin = plot(lambdas*1e9, eta_pL, 'r-o', 'MarkerFaceColor', 'r');
%hold on
%datCoupPSca = scatter(lambdas*1e9, eta_pLs(1,:));
%hold off
ylim(axYlimEtaP)
ylabel('\eta_p [%]', 'FontWeight', fontblAx, 'FontSize', fontszAx)

drawnow

if isSaveGif
    % Save initial frame
    frame = getframe(gifFig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    imwrite(imind,cm,[svfld svnm],'gif', 'Loopcount',inf);
end
else 
    %% Update figure

%-- Update WFE/R 
% Convert WFE from rads to nm
wfr = phz/2/pi*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
% Set out-of-pupil values to NaN for plotting
wfr(wfr==0) = nan;
set(datWFR, 'CData', wfr);

%-- Update Donut PSF
set(datPSF, 'CData', iPSFv_BB);

%-- Update eta_s
    % Note, keep fiber at same location to not "correct" for T/T
%eta_sLs(i,:) = squeeze(eta_mapss(i, eta_sInd(1), eta_sInd(2),:));
eta_sL = eta_maps(eta_sInd(1), eta_sInd(2), :)*eta_sYSCL;
eta_ss(i) = mean(eta_sL);
set(datNull, 'YData', eta_ss);

%-- Update eta_p
eta_pAvg = nan(numWavelengths, N/2-1);
for ii = 1:numWavelengths
    % Azimuthally average
    tmpAvg = VFN_An_radAverage(eta_maps(:,:,ii), [N/2+1, N/2+1]);
    eta_pAvg(ii,:) = tmpAvg;
end
eta_pL = mean(eta_pAvg(:,eta_pInd-prg:eta_pInd+prg),2)*100;
eta_ps(i) = mean(eta_pL);
set(datPeak, 'YData', eta_ps);

%-- Update eta map
set(datEMap, 'CData', eta_maps(:,:,central_band_index));

%-- Update eta_p and eta_s vs. wavelength
set(datCoupSLin, 'YData', eta_sL);
set(datCoupPLin, 'YData', eta_pL);
drawnow

if isSaveGif
    % Capture the plot as an image 
    frame = getframe(gifFig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    imwrite(imind,cm,[svfld svnm],'gif','WriteMode','append'); 
end
%pause(1)

end


end

%% TO-DO LIST!
%{
1) Do zernike decomp of real WF
2) Find good way to show real decomp in the main Figure or in adjacent figure
3) Consider smoothing/re-interpolating the WF data so it's less course
    a) This could be an optional "filtering" feature
%}


%{
    _____________ATTEMPTS AT SAMPLING WF DATA_____________________
  
    % Extract frame
    wfTMP = wfres(:,:,1);
    % Up-sample matrices before rotation:
        % Want enough pts that rot. interp. via nearest-neighbor looks good.
    wfTMP = interp2(wfTMP, 5, 'nearest');
    wfmskR = interp2(wfmsk, 5, 'nearest');
    % Get center of mask in upscaled units (assuming roughly symmetric)
    [rws, cls] = find(wfmskR);  % find min and max row/col indices
    rwc = (max(rws)-min(rws))/2 + min(rws);
    clc = (max(cls)-min(cls))/2 + min(cls);
    % Center up-scaled matrices on mask location
    % Determine rows and columns to add
    vad = round((size(wfTMP)/2-[rwc clc]));
    rwa = zeros(vad(1),size(wfTMP,2));
    cla = zeros(size(wfTMP,1), vad(2));
    if vad(1) > 0
        wfTMP = [rwa; wfTMP];
        wfmskR = [rwa; wfmskR];
    else
        wfTMP = [wfTMP; rwa];
        wfmskR = [wfmskR; rwa];
    end
    if vad(2) > 0
        wfTMP = [cla wfTMP];
        wfmskR = [cla wfmskR];
    else
        wfTMP = [wfTMP cla];
        wfmskR = [wfmskR cla];
    end
    % Rotate centered matrices 
    wfTMP = imrotate(wfTMP, 30);
    wfmskR = imrotate(wfmskR, 30);
    % Get flt2flt on mask for rescaling
        % get furthest point from the center (since already centered)
    [rws, cls] = find(wfmskR);
    f2f = 2*max(sqrt((rws-size(wfmskR,1)/2-1).^2 + (cls-size(wfmskR,2)/2-1).^2));
    % Resample wfTMP at simulated pupil sampling
    wfTMP1 = interp2(0:(2*apRad)/f2f:(size(wfTMP,1)-1)*2*apRad/f2f, ...
                    0:(2*apRad)/f2f:(size(wfTMP,2)-1)*2*apRad/f2f, ...
                    wfTMP, ...
                    0:(size(wfTMP,1)-1)*2*apRad/f2f, ...
                    0:(size(wfTMP,2)-1)*2*apRad/f2f);
                
                
                x = 0:(2*apRad)/f2f:(size(wfTMP,2)-1)*2*apRad/f2f;
y = 0:(2*apRad)/f2f:(size(wfTMP,1)-1)*2*apRad/f2f;
xi = 0:(size(wfTMP,2)-1)*2*apRad/f2f;
yi = 0:(size(wfTMP,1)-1)*2*apRad/f2f;
                [xt yt] = meshgrid(x,y); %create vector arrays
[xit yit] = meshgrid(xi,yi); %create vector arrays
wftmp1 = interp2(xt,yt,wfTMP,xit,yit);
    
    %-- Rotate matrices first and resample to match simulated pupil
    % Apply mask to WF
    wfTMP = wfTMP.*wfmsk;
    % Crop mask and WF to pupil size (for centering purposes)
	[rws, cls] = find(wfmsk);  % find min and max row/col indices
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1:rw2, cl1:cl2);    % crop WF
    wfmskCRP = wfmsk(rw1:rw2, cl1:cl2); % crop mask
    % Rotate both matrices
    wfTMP = imrotate(wfTMP, 30, 'nearest');
    
    % Rotate pupil mask 30deg to match simulated pupil orientation
    wfmskROT = imrotat(wfmsk, 30);
    % Rotate data accordingly
    
    % Extract frame
    wfTMP = wfres(:,:,1);
    % Apply mask
    wfTMP = wfTMP.*wfmsk;
    wfTMP = imrotate(wfTMP, -30, 'nearest');
    % Crop to pupil size
    [rws, cls] = find(imrotate(wfmsk,-30,'nearest'));  % find min and max row/col indices
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1:rw2, cl1:cl2);
    % Resample at simulation pupil diameter
    tmp = interp2(wfTMP, linspace(1,size(wfTMP,2),2*apRad)', linspace(1,size(wfTMP,1),2*apRad), 'nearest');
        
    % Find flat-to-flat value as defined by makeKeckPupil()
        % 1) get radius of pupil
        % 2) get furthest point from the center defined by the radius
    [rws, cls] = find(wfmsk);   % find min and max row/col indices
    rwR = (max(rws)-min(rws))/2;    % radius of pupil along matrix axes
    clR = (max(cls)-min(cls))/2;
    flt2flt  = 2*max(sqrt((rws-rwR-1).^2 + (cls-clR-1).^2))
%}

%{
    ___________ FIRST FUNCTIONAL SAMPLING________________

    % Extract frame
    wfTMP = wfres(:,:,i);
    % Mask using charlotte's pupil
    wfTMP = wfmsk.*wfTMP;
    % Pad sides so that up-sampling produces even-sized grid
    wfTMP = padarray(wfTMP, [1 1]);
    % Up-sample to decrease error sensitivity in centering and rotation
        % Want enough pts that rot. interp. via nearest-neighbor looks good.
    wfTMP = interp2(wfTMP, 5, 'nearest');
    % Recenter and crop (assuming roughly symmetric)
    [rws, cls] = find(wfTMP);  % find min and max row/col indices
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1-1:rw2, cl1-1:cl2);
    % Rotate to fix orientation
    wfTMP = imrotate(wfTMP, 30);
    % Find inscribed circle size
    [rws, cls] = find(~wfTMP);
    crc = 2*min(sqrt((rws-size(wfTMP,1)/2-1).^2 + (cls-size(wfTMP,2)/2-1).^2));
    % OPTIONAL: show inscribed circle
    figure; imagesc(wfTMP)
    viscircles([size(wfTMP,2)/2-1, size(wfTMP,1)/2-1], crc/2);
    % Resize so that inscribed circle is >~= pupCircDiam
    fdg = 5;    % fudge factor on pupCircDiam to make sure bigger than crc
	x = 0:(pupCircDiam+fdg)/crc:(size(wfTMP,2)-1)*(pupCircDiam+fdg)/crc;
    y = 0:(pupCircDiam+fdg)/crc:(size(wfTMP,1)-1)*(pupCircDiam+fdg)/crc;
    xi = 0:(size(wfTMP,2)-1)*(pupCircDiam+fdg)/crc;
    yi = 0:(size(wfTMP,1)-1)*(pupCircDiam+fdg)/crc;
    [xt, yt] = meshgrid(x,y); %create vector arrays
    [xit, yit] = meshgrid(xi,yi); %create vector arrays
    wfTMP = interp2(xt,yt,wfTMP,xit,yit);
    % Pad to match Simulated pupil diameter
    pd = size(wfTMP) - [N N];
    if pd(1) < 0
        % Need to add rows 
        wfTMP = padarray(wfTMP, [-round(pd(1)/2) 0]);
    end
    if pd(2) < 0
        % Need to add columns
        wfTMP = padarray(wfTMP, [0 -round(pd(2)/2)]);
    end
    % Final crop to ensure proper size 
        % Accounts for initial oversize or padding which adds 1 extra row/col
    wfTMP = wfTMP(1:N, 1:N);
    % Apply pupil!
    wfTMP = wfTMP.*PUPIL;
    % OPTIONAL: show final pupil with WF
    figure; imagesc(wfTMP)
%}