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
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model

%--- Define wavefront residuals at the central wavelength
isrealWF = true;      % Flag to mark if real or synthetic WFs should be used
if isrealWF
    % Real WF data parameters
    wfPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
    wfRSN = 'residualWF.fits';          % Filename for residuals
    wfMKN = 'DMPupil.fits';             % Filename for pupil mask on residuals
    % NOTE::: Raw data assumed to be in [Volts]
    % Define sample start and end values
    wfstrt = 100;
    wfSamps = 20;
    % Scaling factor for data (nm/V)
    wfSCL = 600;     
    % Zernikes to plot in figure
    nolls = [2, 3, 5, 6];
    % Flag to choose if WF reconstruction or raw WF should be used
    israwWF = false;
    % Number of zernike modes used in reconstruction
    recNMds = 36;
    %%%% NOTE:::: DEAL WITH MANUAL CORRECTIONS IN WF SECTION BELOW!!!
else
    % Synthetic wavefront parameters
    nolls = [2, 5, 6];
    coeffs = 0*[0.01, 0.025, 0.01];    % Waves RMS at given nolls
    wfSamps = 10;                      % Number of WF samples to synthesize
    isWFPhzMod = true;                % Flag to chose if phase should be modulated
    wfModPhzEnd = 2*pi;                % Temporal modulation phase end (without delay)
end

%-- Define Tip/Tilt residuals  (ADC PARAMS ARE DEFINED SEPARATELY)
% NOTE::: To disable TT residuals, use isrealTT = false and set ttres = [0 0];
isrealTT = true;      % Flag to mark if real or synthetic TTs should be used
% Scaling factor for data (arcsec to waves RMS)
keckD = 10.949;                 % Real-world keck pupil diameter [m] - 10.949 = circumscribed
ttSCL = keckD/206265/lambda0;       % (D in [m])/(arcsec/rad)/(lambda in [m])  = arcsec2wavesPV
% Vectorize to include PV2RMS conversion for both pupil axes separately
    % PV2RMS values determined empirically using generateZernike_fromlist then
    % confirmed by setting ttres = 41.4mas and confirming 1lam/D PSF shift since
    % 41.4mas = 1lam/D @2.2um
ttSCL = ttSCL*ones(2,1)./[3.819513, 4.212690]'; 
if isrealTT
    % Real TT data parameters
    ttPTH = 'C:\Users\danie\OneDrive - California Institute of Technology\Mawet201718\VortexFiberNuller_VFN\Presentations\SPIE_2019\telemetry_20190420\';
    ttRSN = 'residualTT.fits';          % Filename for residuals
    % NOTE::: raw data assumed to be in [arcsec]
    % Define sample start and end values
    if exist('wfstrt', 'var')
        % By default, use wfstrt when defined
        ttstrt = wfstrt;
    else
        % NOTE::: If wfstrt is not defined, user must define a ttsrt manually
        ttstrt = 100;
    end
    % no ttSamps since should sample exactly as many TTs as WFs  
else
    % Synthetic TT parameters
    ttres = 0*[0, 41.4];    % Peak TT error in [mas] at central wavelength
    isTTPhzMod = false;     % Flag to chose if phase should be modulated
    ttModPhzEnd = 2*pi;     % Temporal modulation phase end (without delay) 
end

%-- Define ADC residuals:   in [mas]
% NOTE::: To have disable ADC residuals, set adcres = to all 0's
%adcres = 0*ones(2,numWavelengths);         % Disable ADC residuals
% NOTE::: To define constant ADC residuals, provide 2D matrix here
adcres = 0*[-37.3 -39.4 41.4 43.5 45.6;...         % Tip residuals at each wavelength
            0 0 0 0 0];          % Tilt residuals at each wavelenght
% Make into 3D matrix with dimensionaly: [wfSamps, tip/tilt, wavelength]
adcres = permute(repmat(adcres, 1, 1, wfSamps),[3 1 2]);
% NOTE::: To provide time-varyin ADC:
% Optionally, provide 3D matrix directly with time-varying ADC residuals
    % Dimensionality still needs to match [wfSamps, tip/tilt, wavelength]

% Give offsets for the vortex mask
offsetX = 0*apRad;%0.0952*apRad;
offsetY = 0*apRad;%0.0524*apRad; 

% Flag to plot coupling maps in log scale
islogcoup = true;

% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = false;

%-- Saving parameters
% Flag to save gif
isSaveGif = true;
% Flag to save fits
isSaveFit = true;
% Save folder
svfld = 'C:\Users\danie\Desktop\RealWF_noTT_noADC_Samp20-1020_5Lam\';
% Save name for gif
svnmGif = 'SPIEFig.gif';
% Save name for fits (prefix only, suffix added in code)
svnmFit = 'Dat';

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
    % Fill in closest pixels (causing simulated mask to be severely undersized)
    wfmsk(15,17) = 1; wfmsk(17, 16) = 1; wfmsk(4, 15) = 1;
    
    %-- Aditional nobs to turn for mask features:
    % rot: angle [degrees] by which to rotate data to match our pupil
    % sclfdg: Fudge factor for mask resizing to match simulated pupil to wfmsk
    % cntfdg: Fudge factor for mask centering to match sim pup to wfmsk
    
    %-- Pre-process WF before loop
    % Mask WF using charlotte's pupil
    wfmsk = repmat(wfmsk, [1 1 wfSamps]);
    wfres = wfres.*wfmsk;
    % Pad sides of WF so that up-sampling produces even-sized grid
    wfres = padarray(wfres, [1 1 0]);
    
else
    %*** Synthetic data
    % Create matrix-version of coeffs to simulate changing WFE
    wfMod = ones(wfSamps, length(nolls));
    if isWFPhzMod
        for i = 1:length(nolls)
            % Create random phase delay
            phzShift = rand*2*pi;
            wfMod(:, i) = sin(linspace(0+phzShift, wfModPhzEnd+phzShift, wfSamps))';    % smoothly modulate WFE amplitude
        end
    end
    %wfMod = repmat(wfMod, 1, length(coeffs));
    coeffs = repmat(coeffs, wfSamps, 1);
    coeffs = coeffs.*wfMod;             % Apply modulation
end

%% Display WF Mask (for real data) or Zernikes (for synthetic data)
if isrealWF
    if isPlotSimp
        % Show raw mask
        figure(101);
        imagesc(wfmskRAW);
        axis image;
        title('Original WF Mask - from PyWFS')
        % Show corrected mask
        figure(102);
        imagesc(wfmsk(:,:,1))
        axis image;
        title('User-Corrected WF Mask - from PyWFS')
    end
    
    % NOTE: Zernikes for real data are calculated during main loop
    
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
    title('Zernike Decomposition')
    xlabel('WF Sample')
    ylabel('RMS [waves at \lambda_0]')
    legH = legend(legs);
    title(legH, 'Noll Index')
end    

%% Deal with TT extraction (real data) or synthesis (synthetic data)
% Resize ttSCL to proper dimensionality for ttres
ttSCL = repmat(ttSCL', wfSamps, 1);
if isrealTT
    %*** Real data
    
    %-- Read data 
    ttres = fitsread([ttPTH, ttRSN]);   % Full residuals file
    % Extract specific samples
    ttres = ttres(ttstrt:ttstrt+wfSamps-1,:);
    % Rescale the WFR to waves rms at lambda0
    ttres = ttres.*ttSCL;    
       
else
    %*** Synthetic data
    
    %-- Create matrix-version of coeffs to simulate changing TT
    ttMod = ones(wfSamps, length(ttres));
    if isTTPhzMod 
        % Modulate phase as requested by user
        for i = 1:length(ttres)
            % Create random phase delay
            phzShift = rand*2*pi;
            ttMod(:, i) = sin(linspace(0+phzShift, ttModPhzEnd+phzShift, wfSamps))';    % smoothly modulate TT amplitude
        end
    end
    ttres = repmat(ttres, wfSamps, 1);
    ttres = ttres.*ttMod;             % Apply modulation
    %-- Convert from mas to [waves RMS]
    ttres = ttres*1e-3.*ttSCL;        % *1e-3 to go from mas to arcsec for ttSCL
end

%-- Final formatting 
% Convert to 3D for adc combination. New dim.: [wfSamps, tip/tilt, wavelength];
ttres = repmat(ttres, 1, 1, numWavelengths);
% Scale tilt by wavelength
ttlamscl = permute(repmat(lambdas/lambda0,wfSamps,1,2),[1 3 2]);
ttres = ttlamscl.*ttres;


%% Display TT residuals (no atmospheric dispersion yet)
ttFig = figure(104);
% Set axes colors to black [0 0 0]
set(ttFig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
yyaxis right
% Plot right axis [waves RMS]
ttPLT = plot(1:wfSamps, ttres(:,1,1), 'r-o', 1:wfSamps, ttres(:,2,1), 'b-o', 'MarkerSize', 4, 'LineWidth', 2);
title('Tip Tilt Residuals - Before Atm. Disp.')
xlabel('WF Sample')
legend('Tip', 'Tilt');
ylabel('RMS [waves at \lambda_0]')
% Save ylimits for rescaling on left side
ttYlimL = ylim;
% Add left label [mas]
yyaxis left
ylabel('RMS [mas]')
% Rescale using mean of ttSCL to account for difference in Tip/Tilt PV2RMS
ylim(ttYlimL*1e3/mean(ttSCL(1,:)))

%% Deal with atmospheric dispersion
% Resize ttSCL to proper dimensionality for ttres
ttSCL = repmat(ttSCL, 1, 1, numWavelengths);

% Convert from [mas] to [waves]
adcres = adcres*1e-3.*ttSCL;

% Add in ADC residuals as additional TT errors directly
ttres = ttres+adcres;

%% Make vortex mask (once before loop since does not change)
EPM = generateVortexMask( charge, coords, [offsetX offsetY] );

%% Generate fibermode at each lambda (once before loop since does not change)
    % Generating fib modes once improves runtime by 25% for N=2048 and wfSamps=10
    
% Parameters for Thorlabs SM600
    % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 11e-6/2;% Core radius [um]
fiber_props.n_core = 1.4606;% core index (interpolated from linear fit to 3 points)
fiber_props.n_clad = 1.4571;% cladding index (interpolated from linear fit to 3 points)
Fnum = 5; % focal ratio of the beam at the fiber
fiber_props.type = 'bessel';

% Iterate through wavelengths generating modes
fibmodes = nan(N, N, numWavelengths);
for ch = 1:numWavelengths
    fibmodes(:,:,ch) = generateSMFmode(fiber_props, lambdas(ch), Fnum, lambdaOverD, coords);
end

%% Generate unit zernike modes (once before loop since does not change)
    % Generating unit modes once improves runtime by 25% on top of fiber mode
    % improvement for N=2048 and wfSamps=10

% Define zernikes on which to decompose
if israwWF
    %-- Reconstruction not requested; use raw WF
    % Set loop to only decompose requested modes (nolls)
    jvls = nolls;
else
    %-- Use reconstruction for WF instead of raw values
    % Set loop to decompose on all recNMds
    jvls = 1:recNMds;
end

% Preallocate zernike matrix
unitzrn = nan(N, N, length(jvls));
% Create basis
for j = jvls
    % Generate zernike (on keck pupil) with 1 wave RMS error
    unitzrn(:,:,j) = generateZernike_fromList(j, 1, PUPIL, pupCircDiam/2, coords);
    % Rescale zernike from radians to waves
    unitzrn(:,:,j) = unitzrn(:,:,j)/(2*pi);
end

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
axYlimEtaS = [0 2];      % y axis limit (in units of eta_sYSCL)
eta_sYSCL  = 1e3;           % scaling for y axis
% realWF Zernike plot
zrnSCL = 1e2;               % scaling factor for y axis
axYlimZrn = [-10 10];          % y axis limits (in units of 1e3)

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

if isrealWF
    %-- Preallocate realWF zernike coeffs FOR PLOTTING
        % NOTE: WF reconstruction is still done over 36 zernikes
    coeffs = nan(wfSamps, length(nolls));
end

% Determine central band 
central_band_index = ceil(numWavelengths/2);

for i = 1:wfSamps
fprintf('wfSamp %03d of %03d\n', i, wfSamps);
    
%% Define pupil field
%*** DEAL WITH WAVEFRONT RESIDUALS
if isrealWF  
    
    %-- FORMAT REAL WAVEFRONT FOR SIMULATION
    % Extract frame
    wfTMP = wfres(:,:,i);
    % Up-sample to decrease error sensitivity in centering and rotation
        % Want enough pts that rot. interp. via nearest-neighbor looks good.
    wfTMP = interp2(wfTMP, 5, 'nearest');
    % Recenter and crop (assuming roughly symmetric)
    [rws, cls] = find(wfTMP);  % find min and max row/col indices
    cntfdg = [0 -10];   % centering fudge factor  [row col] ** must be even
    rw1 = min(rws); rw2 = max(rws);
    cl1 = min(cls); cl2 = max(cls);
    wfTMP = wfTMP(rw1-1+cntfdg(1)/2:rw2+cntfdg(1)/2, cl1-1+cntfdg(2)/2:cl2+cntfdg(2)/2);
    % Rotate to fix orientation
    rot = 27;       % 27 provides decent alignment for furthest segments
    wfTMP = imrotate(wfTMP, rot);
    % Find inscribed circle size
    [rws, cls] = find(~wfTMP);
    crc = 2*min(sqrt((rws-size(wfTMP,1)/2-1).^2 + (cls-size(wfTMP,2)/2-1).^2));
    % OPTIONAL: show inscribed circle
    %figure; imagesc(wfTMP)
    %viscircles([size(wfTMP,2)/2-1, size(wfTMP,1)/2-1], crc/2);
    %title('Inscribed Circ B4 Resize')
    % Resize so that simulated pupil closely matches real pupil (wfmsk)
    sclfdg = -23;    % fudge factor for pupil resizing ** (-)=grow mask; (+)=shrink mask
	x = 0:(pupCircDiam+sclfdg)/crc:(size(wfTMP,2)-1)*(pupCircDiam+sclfdg)/crc;
    y = 0:(pupCircDiam+sclfdg)/crc:(size(wfTMP,1)-1)*(pupCircDiam+sclfdg)/crc;
    xi = 0:(size(wfTMP,2)-1)*(pupCircDiam+sclfdg)/crc;
    yi = 0:(size(wfTMP,1)-1)*(pupCircDiam+sclfdg)/crc;
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
    % Apply pupil!  (Not actually, pupil will be applied at Epup calc)
    % wfTMP = wfTMP.*PUPIL;
    % OPTIONAL: show final pupil with WF
    %figure; imagesc(wfTMP.*PUPIL); title('Final Raw WF - Simulated Pupil')
    
    %-- DECOMPOSE WAVEFRONT INTO ZERNIKES;
    recCof = nan(recNMds,1);
    if ~israwWF
        % Preallocate reconstructed WF matrix
        recWF = zeros(N);
    end
    for j = jvls
        % Extract current mode
        tmpzrn = unitzrn(:,:,j);
        % Left-divide to get coefficient
        recCof(j) = tmpzrn(logical(PUPIL))\wfTMP(logical(PUPIL));
        if ~israwWF
            recWF = recWF + tmpzrn*recCof(j);
        end
    end
    coeffs(i,:) = recCof(nolls);
    
    %-- DISPLAY REQUESTED ZERNIKE COEFFICIENTS (nolls)
    if i == 1
        % Show requested coefficients
        zrnFig = figure(100);
        legs = cell(length(nolls),1);
        for j = 1:length(nolls)
            cofPlt(j) = plot(coeffs(:,j)*zrnSCL, '-o', 'MarkerSize', 4, 'LineWidth', 2);
            if j == 1
                hold on
            end
            legs(j) = {num2str(nolls(j))};
        end
        hold off
        xlim([0 wfSamps]);
        ylim(axYlimZrn);
        title('Zernike Decomposition (Approx in Pupil)', 'FontSize', fontszTi);
        xlabel('WF Sample')
        ylabel('RMS [% waves at \lambda_0]')% \times 10^{-2}]')
        legH = legend(legs);
        title(legH, 'Noll Index')
    else
        for j = 1:length(nolls)
            set(cofPlt(j), 'YData', coeffs(:,j)*zrnSCL);
        end
    end
    
    %-- SET WF FOR USE AND CONVERT FROM WAVES TO RADIANS
    if israwWF
        wfphz = wfTMP*2*pi;
    else
        % Renormalize one last time to ensure RMS values match.
        rmsWF = sqrt(mean(wfTMP(logical(PUPIL)).^2));
        recWF = rmsWF*recWF/sqrt(mean(recWF(logical(PUPIL)).^2));
        wfphz = recWF*2*pi;
    end
    
    %-- DISPLAY PUPIL PROJECTIONS
    if isPlotSimp && (i == 1)
        crcPup = makeCircularPupil(pupCircDiam/2, N);
        wfPLT = wfTMP; wfPLT(wfPLT==0) = nan;
        figure(103); 
        imagesc(xvals/(pupCircDiam/2), yvals/(pupCircDiam/2), wfPLT+0.05*PUPIL+0.1*crcPup); 
        axis image
        axis([-1.2 1.2 -1.2 1.2])
        title('Pupil Projections')
    end
    
else
    %-- Apply synthetic wavefronts
    wfphz = generateZernike_fromList( nolls, coeffs(i,:), PUPIL, pupCircDiam/2, coords); 
end

if isPlotSimp
    figure(2);
end
for ch = 1:numWavelengths
    %*** DEAL WITH TT AND ADC RESIDUALS (in loop to do by wavelength)
    % Check if there are any non-zero values in ttres in case user has disabled 
        % TT/ADC residuals via setting all to 0
    if find(ttres(:,:,ch))
        %-- Non-zero value found thus apply TT/ADC residuals
        % Build TT/ADC wavefronts
        ttphz = generateZernike_fromList([2 3], ttres(i,1:2,ch), PUPIL, pupCircDiam/2, coords);
    else
        ttphz = zeros(N, N);
    end
    
    %*** COMBINE ALL PHASE ERRORS 
    pupphz = ttphz + lambda0/lambdas(ch)*wfphz;
    Epup(:,:,ch) = exp(1i*pupphz).*PUPIL;
    
    if isPlotSimp
        subplot(1,numWavelengths,ch);
        imagesc(xvals/(pupCircDiam/2),yvals/(pupCircDiam/2),angle(Epup(:,:,ch)));
        axis image; 
        axis([-1 1 -1 1]);
        title(['Phase at ',num2str(lambdas(ch)*1e9),'nm']);
        colorbar; 
        colormap(parula(256));
    end
drawnow;   
end

%% Check PSF at each wavelength
% Calculate PSF at each wavelength and display. If titls are set to 1lam/D, PSF
% at each wavelength should also be 1lam/D
% myfft2 = @(x) fftshift(fft2(fftshift(x)));
%     
% iPSF = zeros(coords.N,coords.N,length(lambdas)); 
% 
% ch = 1;% channel index
% for lam = lambdas 
%     PSF = myfft2(Epup(:,:,ch)); % PSF (complex field)
%     iPSF_lam = abs(PSF).^2/normI;
%     lam_frac = lam/lambda0;
%     iPSF(:,:,ch) = (1/lam_frac)^2*interp2(coords.X,coords.Y,iPSF_lam,coords.X/lam_frac,coords.Y/lam_frac,'linear',0);
%     ch = ch + 1;
% end
% 
% figure(9);
% for ch = 1:numWavelengths
%     subplot(1,numWavelengths,ch);
%     if islogcoup
%         imagesc(xvals/(lambdaOverD*lambdas(ch)/lambda0),yvals/(lambdaOverD*lambdas(ch)/lambda0),log10(iPSF(:,:,ch)));
%     else
%         imagesc(xvals/(lambdaOverD*lambdas(ch)/lambda0),yvals/(lambdaOverD*lambdas(ch)/lambda0),iPSF(:,:,ch));
%     end
%     axis image;
%     axis([-2 2 -2 2]);
%     title(['\eta at ',num2str(lambdas(ch)*1e9),'nm']);
%     colorbar;
%     if islogcoup; caxis([-4 0]); end
%     colormap(gray(256));
% end

%% Get PSF without vortex mask

if isPlotSimp
    % Get broadband PSF; used only for this plot so not necessary
    iPSF_BB = getPSF(Epup,lambda0,lambdas,normI,coords);

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
if isPlotSimp
    figure(4)
    imagesc(xvals/(pupCircDiam/2),yvals/(pupCircDiam/2),angle(EPM(:,:,central_band_index).*Epup(:,:,central_band_index)));
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
% Use pre-made fibermodes to improve runtime
eta_maps = generateCouplingMap_polychromatic(Epup.*EPM, fiber_props, lambda0, Fnum, lambdas, totalPower0, lambdaOverD, 3*lambdaOverD, coords, fibmodes);

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
wfr = wfphz/2/pi*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
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
if isSaveFit
    fitswrite(wfr, [svfld sprintf('wfr%06d.fits',i)]);
end

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
if isSaveFit
    fitswrite(iPSFv_BB, [svfld sprintf('psfBB%06d.fits',i)]);
end

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
if isSaveFit
    fitswrite(eta_maps, [svfld sprintf('etaMaps%06d.fits',i)]);
end

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
    imwrite(imind,cm,[svfld svnmGif],'gif', 'Loopcount',inf);
end
else 
    %% Update figure

%-- Update WFE/R 
% Convert WFE from rads to nm
wfr = wfphz/2/pi*lambda0*1e9.*PUPIL;      %.*PUPIL Optional to show with pupil
% Set out-of-pupil values to NaN for plotting
wfr(wfr==0) = nan;
set(datWFR, 'CData', wfr);
if isSaveFit
    fitswrite(wfr, [svfld sprintf('wfr%06d.fits',i)]);
end

%-- Update Donut PSF
set(datPSF, 'CData', iPSFv_BB);
if isSaveFit
    fitswrite(iPSFv_BB, [svfld sprintf('psfBB%06d.fits',i)]);
end

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
if isSaveFit
    fitswrite(eta_maps, [svfld sprintf('etaMaps%06d.fits',i)]);
end

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
    imwrite(imind,cm,[svfld svnmGif],'gif','WriteMode','append'); 
end
%pause(1)

end


end

%% TO-DO LIST!
%{
1) Do zernike decomp of real WF
    *** DONE
2) Find good way to show real decomp in the main Figure or in adjacent figure
3) Consider smoothing/re-interpolating the WF data so it's less course
    *** DONE - reconstructing WF using zernike decomp
    a) This could be an optional "filtering" feature
4) Find a way to fit more of our simulated pupil within the real pupil
    Options:
    *** DONE
    a) Rotate WF less (<30) so that orientation is more favorable
    b) Don't use wfmsk at all for sizing, only for centering
5) Display mods to original pupil as different color in figure 101
6) Add slide of raw DM voltages without mask provided by charlotte
    *** DONE
    - Say that this is where we get our WF data from
7) Prevent simplePlot figures from being recast each time
    - Modify so that they are only generated once and then only the data gets
        updated on each iteration in the same way that the SPIE fig is handled.
8) Consider scaling zernike coefficients by rmsWF/rmsRec ratio
    - This may be a valid way of rescaling the coefficients to get more accurate
        reported values.
9) Add ADC errors
    *** DONE
    - separate section from TT
    - Net TTcoeffs matrix should be 3D [sample, dir(T/T), wavelength]
10) T/T residuals
    *** DONE
11) Calc planet coupling at star location, not center of frame
    - w/ TT and ADC residuals, the star is no longer at the center of the frame
    - Thus, we should calc. the lam/D offset due to TT and ADC and then center
    the planet averaging around this point.
12) Modify generateCouplingmap_polychromatic.m to improve efficiency
    *** DONE
    - 25% of runtime is currently spent in this function
    - Specifically, most of the time is spent generating the bessel SMF model
    - This only needs to be generated once, not each time. As such, modify code
        to generate bessel model once and then re-use it each time.
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