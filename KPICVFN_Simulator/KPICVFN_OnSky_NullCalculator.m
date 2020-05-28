%{
Code for calculating the null depth of a VFN observation given a set of WFE.

*** Based heavily off KPICVFN_OnSky_EMapGenerator.m

*** This script is meant to be used with KPICVFN_OnSky_WFDecomposer

Goals:___
- Take in WF residual maps and determine the corresponding null depth
- Optionally take in a list of zernikes to create WF residuals
- Allow user to vary the amount of each zernike used in creating WF

Notes:___
* Does not do anything with T/T or ADC residuals at the moment.
* Does not calculate planet coupling; only deals with null point
* The final RMS within the PUPIL will not be equal to the RSS of the zernikes
    - in theory you can RSS the zernike coeffs to get the total RMS
    - this does not hold on non-circular pupils though since zernikes are only
        orthogonal on the unit circle
    ** As such, I save the final RMS on the PUPIL for future reference
* Does not renormalize the wfr to set RMS of wfr to the rms from the WF decomp
    - since the user will be modulating the zernikes individually, this
        renormalization doesn't make much sense in most cases
    ** If you want the reconstructed wavefront that is most equivalent to
        Charlotte's data (ie renormalzied or even the raw one) use either
        wfr_reform or wfr_recons files and set isZernWFR false. This will allow
        you to simulate the actual KPIC performance directly.
%}

clear; close all; 
addpath('C:\Users\danie\Documents\MATLAB\VFN-Simulations\VFNlib');
addpath(genpath('C:\Users\danie\Documents\MATLAB\VFN-Lab\AnalysisCode\'))
addpath(genpath('C:\Users\danie\Documents\MATLAB\falco-matlab'))

%% Input parameters 

%-- Extract values below from fits header
% Folder with files
wfrfld = 'C:\Users\danie\Documents\VFN_Simulations_Temporary\FirstBigRun\';
% filename counter format
wfrFMT = '%06d';
% fits to reference
fitsnm = sprintf(['wfr_reform' wfrFMT '.fits'],1);
% Extract keywords
kwds    = fitsinfo([wfrfld fitsnm]);
kwds    = kwds.PrimaryData.Keywords;      % all keywords


%-- Provide regular parameters
% Define smapling info
N = VFN_An_getKwd(kwds, 'N'); % Size of computational grid (NxN samples) 
apRad = VFN_An_getKwd(kwds, 'aprad'); % Aperture radius in samples 

% Define wavelength info
lambda0 = VFN_An_getKwd(kwds, 'lambda0'); %central wavelength
fracBW = 0.20; %\Delta\lambda/\lambda
numWavelengths = 5;% number of discrete wavelengths 
lambdas = getWavelengthVec(lambda0,fracBW,numWavelengths);% array of wavelengths (meters)

%Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
f = 25e-3;      %[m] final focal length
D = 13e-3;      %[m] pupil size at DM
    % Since using falco MFT propogator, can choose our final image sampling
lambdaOverD_samp = 10; %[samples/ (lam/D)] in units of samples in the image plane
im_size = 6;    %[lam0/D] field of view in final image plane

% Define charge of the vortex mask at each wavelength
charge = 2*ones(1,numWavelengths); % achromatic
%charge = 2*lambda0./lambdas; % simple scalar model


%-- Define WFR 
% WFR filename (can be fits of frames OR fits of zernike coeffs)
    % Use isZernWFR to denote which one it is
wfrfnm = 'zern_coeff'; %wfr_recons
rmsfnm = 'rms_WF';      % filename for rms values (col1 = on charlotte pupil; col2 = on Keck PUPIL)
% Flag to use Zern Cofs to generate wfr. true=generate from zernikes, false=use pregenerated by decompose code
isZernWFR = true;
% Pupil filename (desired pupil should be 3rd frame of fits cube)
pupfnm = 'wfr_pupils';
%-- Define modifications to WFR 
% Factor by which to modulate the given zernike mode (should be fractional)
    % Must be (1,N) dimensionality, NOT (N,1)
zrnMod = ones(1,VFN_An_getKwd(kwds, 'RECNMDS'));
zrnMod(5:6) = 1;        % Sample increase of 2x primary astig
zrnMod(2:3) = 0;        % remove T/T from this sim since it'll change PSF centering
% *** Read note header about how to get the most accurate WFR reconstruction
% with the correct renormalization (for when trying to simulate the real KPIC)

%-- Get total number of samples from keywords
wfSamps = VFN_An_getKwd(kwds, 'wfsamps');  % Old method: size(dir([wfrfld wfrfnm '*.fits']),1);


%--Offsets for the vortex mask
offsetX = 0*apRad;
offsetY = 0*apRad;


%-- Control Flags
% Flag to use fits file pupil. true=fits, false=regenerate and rotate pupil
isFitsPup = true;
% Flag to plot intermediate figures (before SPIE figure)
isPlotSimp = true;
% Flag to make figures visible
visflg = 'off';     % should be 'on' or 'off'
% Flag to make all figures visible at the end
isFigsVisEnd = true;


%-- Saving parameters
% Flag to save gif
isSaveGif = true;
% Flag to save fits
isSaveFit = true;
% Save folder
svfld = [wfrfld 'ZernNoMod_Charge' num2str(charge(1)) '/'];
if isSaveFit
    % Create foler if it doesn't already exist
    mkdir(svfld)
end
% Save name for gif
svnmGif = 'PSF_EMaps.gif';
% Delay time for gif (time between frames in [s])
gifDelay = 0.1;
% Checkpoints for save data
    % Number of time instances to save per fits file
SAVE_CHK = VFN_An_getKwd(kwds, 'naxis3');

%% Generate the coordinate system

%-- Coordinates in the focal plane
coords = generateCoordinates(N);% Creates NxN arrays with coordinates 
xvals = coords.xvals;% Helpful for plotting
yvals = coords.yvals;

%-- Coordinates in the pupil plane
Nxi = im_size*lambdaOverD_samp;
coordsfp = generateCoordinates( Nxi );

%-- Key values for getting scaling right using falco MFT propagator
% Lambda over D of central wavelength in radians at final focal plane 
lambdaOverD_rad = lambda0*f/D;     
% "pixel" size in pupil and image planes
dx  = D/(2*apRad);
dxi = lambdaOverD_rad/lambdaOverD_samp;


%% Create array with pupil function

if isFitsPup
    PUPIL = fitsread([wfrfld pupfnm '.fits']);
    PUPIL = PUPIL(:,:,3);
    % Assume corresponding wfr_reform was provided for keywords in header
    pupCircDiam = VFN_An_getKwd(kwds, 'pupCrcDm');  
    %-- Pad back to NxN
    PUPIL = pad_crop(PUPIL,N);
else
    [PUPIL, pupCircDiam] = makeKeckPupil( 2*apRad, N );

    % Rotate pupil to match orientation of PyWFS data
    % Assume corresponding wfr_reform was provided for keywrds in header
    PUP_rot = VFN_An_getKwd(kwds, 'pup_rot');
    PUPIL = imrotate(PUPIL, PUP_rot, 'crop');    % 'crop' maintains size and centering
end

%-- Get normalization factors
% Norm for coupling fractions (simple sum since using MFT propagator)
totalPower0 = sum(abs(PUPIL(:)));
% Norm for PSF plots (peak of ideal PSF)
PSF = propcustom_mft_PtoF(PUPIL, f, lambda0, dx, dxi, Nxi, dxi, Nxi);
normI = max(abs(PSF(:)).^2);


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
    
% % Parameters for Thorlabs SM600
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
% fiber_props.core_rad = 11e-6/2;% Core radius [um]
% fiber_props.n_core = 1.4436;% core index 
% fiber_props.n_clad = 1.4381;% cladding index 
% Fnum = getMFD(fiber_props,lambda0)/(lambda0*1.4); % focal ratio of the beam at the fiber
% fiber_props.type = 'bessel';

% Iterate through wavelengths generating modes
fibmodes = nan(Nxi, Nxi, numWavelengths);
for ch = 1:numWavelengths
    fibmodes(:,:,ch) = generateSMFmode_gaussian(1.4*(lambdas(ch)/lambda0)*lambdaOverD_samp,coordsfp);
    % Old method:
    %fibmodes(:,:,ch) = generateSMFmode(fiber_props, lambdas(ch), Fnum, lambdaOverD, coords);
end

%% Load wfr map info if not generating them from zernikes each time
if ~isZernWFR
    %-- Get keywords which describe shape of wfr files
    kwdsN    = fitsinfo([wfrfld sprintf([wfrfnm wfrFMT '.fits'],1)]);
    kwdsN    = kwdsN.PrimaryData.Keywords;      % all keywords
    frames_per_file = VFN_An_getKwd(kwdsN, 'NAXIS3');
    
%   %-- Following block of code isn't needed since wfr is read in dynamically in
%       % main for loop
%     % Preallocate full file
%     wfr = nan(VFN_An_getKwd(kwdsN, 'NAXIS1'), VFN_An_getKwd(kwdsN, 'NAXIS2'), wfSamps);    
%     
%     %-- Read all files and place in wfr
%     % Get total number of files over which to iterate
%     N_wfr = size(dir([wfrfld wfrfnm '*.fits']),1);
%     for i = 1:N_wfr
%         % read the file
%         wfrTMP = fitsread([wfrfld sprintf([wfrfnm wfrFMT '.fits'],i)]);
%         % Place in matrix
%         if i == N_wfr
%             % On the last file so change indexing since it may have less frames
%             wfr(:,:,(frames_per_file*(i-1))+1:end) = wfrTMP;
%         else
%             % On every other file, simply place it in the cube
%             wfr(:,:, (frames_per_file*(i-1))+1:frames_per_file*i) = wfrTMP;
%         end
%     end  
% 
%     %-- Check if any frames were left unallocated
%     % find the number of nan values in each frame
%     nan_vec = squeeze(sum(isnan(wfr),[1 2]));
%     % check for missed frames
%     if sum(nan_vec(nan_vec == size(wfr,1)*size(wfr,2))) ~= 0
%         error('A wfr frame was left unallocated')
%     end
% 
%     %-- Change nan's back to 0's
%     wfr(isnan(wfr)) = 0;
%     
%     %-- Reclaim some memory
%     clear wfrTMP kwdsN 
end


%% Prepare for using zernikes to create wfr
if isZernWFR
    %-- Preallocate zernike and rms matrices
    % Get keywords which describe shape of zernike files
    kwdsN    = fitsinfo([wfrfld sprintf([wfrfnm wfrFMT '.fits'],1)]);
    kwdsN    = kwdsN.PrimaryData.Keywords;      % all keywords
    frames_per_file = VFN_An_getKwd(kwdsN, 'NAXIS2');
    % Preallocate full file
    zrn_cof = nan(wfSamps, VFN_An_getKwd(kwdsN, 'NAXIS1')); 
    rmsWF = nan(wfSamps, 2);        % use 2 since we know rms matrix has 2 cols
    
    %-- Read all files and place in zrnike/rmsWF matrix
    % Get total number of files over which to iterate
    N_zrn = size(dir([wfrfld wfrfnm '*.fits']),1);
    for i = 1:N_zrn
        % read the files
        zrn_tmp = fitsread([wfrfld sprintf([wfrfnm wfrFMT '.fits'],i)]);
        rms_tmp = fitsread([wfrfld sprintf([rmsfnm wfrFMT '.fits'],i)]);
        % Place in matrix
        if i == N_zrn
            % On the last file so change indexing since it may have less frames
            zrn_cof((frames_per_file*(i-1))+1:end, :) = zrn_tmp;
            rmsWF((frames_per_file*(i-1))+1:end, :) = rms_tmp;
        else
            % On every other file, simply place it in the matrix
            zrn_cof((frames_per_file*(i-1))+1:frames_per_file*i, :) = zrn_tmp;
            rmsWF((frames_per_file*(i-1))+1:frames_per_file*i, :) = rms_tmp;
        end
    end  
    
%     %-- Check if any frames were left unallocated
%     if any(isnan(zrn_cof),'all')
%         error('A zernike coeff value was left unallocated')
%     end

    %-- Modulate coeffs by user-requested values
    zrn_cof = zrnMod.*zrn_cof;

    %-- Create unit zernikes
    % Set loop to create as many modes as are in zrn_cof data
    jvls = 1:size(zrn_cof, 2);

    % Preallocate zernike matrix
    unitzrn = nan(N, N, length(jvls));
    % Create basis
    for j = jvls
        % Generate zernike (on keck pupil) with 1 WAVE RMS error
        unitzrn(:,:,j) = (1/(2*pi))*generateZernike_fromList(j, 1, PUPIL, pupCircDiam/2, coords);
    end

    if any(isnan(unitzrn),'all')
        error('an element in unitzrn was left unallocated')
    end
    
    %-- Save the zrnMod vector for future reference
    if isSaveFit
        fitswrite(zrnMod, [svfld 'zenMod.fits']);
    end
    
    %-- Reclaim some memory
    clear zrn_tmp rms_tmp
end

%% Prepare to loop through WF samples
% Use subtightplot for control over plot spacing
% subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.06], [0.1 0.05], [0.05 0.05]);

%-- Define font and linewidth sizes
fontszAx = 14;
fontblAx = 'normal';
fontszTi = 14;

%-- Axis limits
% PSF plot
axlimPSF = [-3 3 -3 3];     % x-y axis limits
% Eta Map
% axlim2EMap = axlimPSF;      % x-y axis limits
% Eta_p plots
% axYlimEtaP = [0 20];        % y axis limit

%-- Coloraxes
% PSF
    % Determined empirically from several samples
    % Feel free to modify as needed
caxlimPSF = [0 0.33];  
% Eta map
% caxlimEMap = axYlimEtaP/100;


% Determine central band 
% central_band_index = ceil(numWavelengths/2);

%% actual loop
% Preallocate null depth matrix
eta_ss = nan(wfSamps, numWavelengths);
% Preallocate broadband PSF image matrix
iPSF_BB = nan(Nxi, Nxi, SAVE_CHK);
if isZernWFR
    % Preallocate the rms of the final WF after modulation and reconstruction
        % only done when isZernWFR since otherwise the rms is the same as rmsWF
    rmsrecons = nan(wfSamps,1);
end

for i = 1:wfSamps
fprintf('wfSamp %06d of %06d (frame %3d in cube %3d)\n', i, wfSamps, mod(i-1,SAVE_CHK)+1, floor(i/SAVE_CHK)+1 - ~mod(i,SAVE_CHK));
    
%%  Define pupil field

%-- Format WFR map as needed
if isZernWFR
    %-- Generate wfr map from pre-loaded zernikes
    % Use vectorized method to reconstruct (runs 2x as fast as for-loop)
        % reshape zernike vector into size = [1 1 3];
        % Multiply unitzrn by this vector elementwise (matlab implicitly replicates the vector so that dimensions agree)
        % Sum along 3rd access to get final reconstruction
    wfr = sum(unitzrn.*reshape(zrn_cof(i,:),1,1,size(zrn_cof,2)), 3);
    
    
    %-- Calculate rms on reconstructed wfr (including mods from zrn_mod)
    rmsrecons(i) = rms(wfr(logical(PUPIL)));
    
    
%     % Following chunk of code is not done anymore. See header for details.
%     %-- Rescale so that final rms matches original rms before reconstruction
%         % If we're modulating the the zernikes, we probably shouldn't force the rms to match since
%         %       our goal is to see if the net rms is a good proxy for null depth
%     wfr = wfr*(rmsWF(i,2)/sqrt(mean(wfr(logical(PUPIL)).^2)));
    
    %-- Convert from waves to phase
    wfr = wfr*2*pi;
else
    % *** Using saved wfr maps
    %-- Load new cube if necessary
%     % (following bit of commented code verifies the indexing used in the following lines)
%     for i = 1:18
%     fprintf('itr %d, new? %d, cbind %d, frind %d\n',i, ~mod(i-1, 5), floor(i/5)+1 - ~mod(i,5), mod(i-1,5)+1)
%     end
    if ~mod(i-1, frames_per_file)
        % read in new cube of WF data
        wfrCB = fitsread([wfrfld sprintf([wfrfnm wfrFMT '.fits'], floor(i/frames_per_file)+1 - ~mod(i,frames_per_file))]);
        % Change nan's back to 0's
        wfrCB(isnan(wfrCB)) = 0;
    end
    
    %-- Extract frame from existing cube
    wfr = wfrCB(:,:, mod(i-1,frames_per_file)+1);
    
    %-- Pad wfrTMP to full size
    wfr = pad_crop(wfr, N);
    
    %-- Convert wfr from nm to phase
    wfr = wfr*2*pi/lambda0/1e9;
end

%% Calculate PSF
for ch = 1:numWavelengths
    % create electric field at pupil (including vortex)
    Epupv = exp(1i*wfr*lambda0/lambdas(ch)).*PUPIL.*EPM(:,:,ch);   % also rescale to lambda
    
    % get PSF 
    PSFv(:,:,ch) = propcustom_mft_PtoF(Epupv,f, lambdas(ch), dx, dxi, Nxi, dxi, Nxi);
end

%-- Get broadband image
% Reset image matrix whenever we've completed a SAVE_CHK checkpoint
if ~mod(i-1, SAVE_CHK)
    iPSF_BB = iPSF_BB*nan;
end
% Calculate intensity (peak normalized)
iPSF_tmp = abs(PSFv).^2/normI;
% Take mean and store in image matrix
iPSF_BB(:,:,mod(i-1,frames_per_file)+1) = mean(iPSF_tmp,3);

%% Calculate null 
for ch = 1:numWavelengths
    eta_ss(i,ch) = (abs(sum(sum(PSFv(:,:,ch).*fibmodes(:,:,ch)))).^2)/totalPower0;
end

%% Plot if desired
if isPlotSimp
    %-- Display PSF (Broadband) and EMap (lambda0) (top row in SPIE figure)
    if i == 1
        %-- First iteration: make the figure
        fig = figure(7);
        set(fig, 'units', 'inches', 'Position', [0 0 6 6], 'color', 'w', 'Visible', visflg);
        %-- Plot instantaneous Donut PSF (polychromatic)
%         axPSF = subplot(1,2,1);
        datPSF = imagesc(coordsfp.xvals/lambdaOverD_samp,coordsfp.yvals/lambdaOverD_samp,iPSF_BB(:,:,mod(i-1,frames_per_file)+1));
        axis image; 
        set(datPSF.Parent, 'TickDir', 'out');
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
        colormap(datPSF.Parent, gray(256));
        caxis(caxlimPSF)
        
%         %-- Plot instantaneous EMap (at lambda0)
%         axEMap = subplot(1,2,2);
%         datEMap = imagesc(xvals/lambdaOverD,yvals/lambdaOverD,eta_maps(:,:,central_band_index));
%         % Can't show circle where eta_p will be calculated since line-profiles
%         % are not generated in this script, thus no peak is calculated
%         axis image;
%         axis(axlim2EMap);
%         set(axEMap, 'TickDir', 'out');
%         title(['Coupling at ',num2str(lambda0*1e9),'nm'], 'FontSize', fontszTi);
%         xlabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
%         ylabel('\lambda/D', 'FontWeight', fontblAx, 'FontSize', fontszAx)
%         EMapCBAR = colorbar; 
%         EMapCBAR.Label.String = '\eta';
%         EMapCBAR.Label.FontWeight = fontblAx;
%         EMapCBAR.Label.FontSize = fontszAx;
%         %set(axEMap, 'Colormap', gray(256))
%         colormap(axEMap, gray(256));
%         caxis(caxlimEMap)
    else
        %-- Other iterations:
        % Update Donut PSF
        set(datPSF, 'CData', iPSF_BB(:,:,mod(i-1,frames_per_file)+1));
%         set(datEMap, 'CData', eta_maps(:,:,central_band_index));
    end

    drawnow;
end


if isSaveFit && ~(mod(i,SAVE_CHK))
    %-- Save fits files
    fitswrite(eta_ss(i-SAVE_CHK+1:i,:), [svfld sprintf('etaStar%06d.fits',floor(i/SAVE_CHK))]);
    fitswrite(iPSF_BB, [svfld sprintf('psfBB%06d.fits',floor(i/SAVE_CHK))]);
    if isZernWFR
        % save rms of final wfr used after reconstruction with modulation
        fitswrite(rmsrecons(i-SAVE_CHK+1:i), [svfld sprintf('rmsrecons%06d.fits',floor(i/SAVE_CHK))]);
    end
elseif isSaveFit && i==wfSamps
    %-- Save final fits in case it is not an integer multiple of SAVE_CHK
    fitswrite(eta_ss(i-mod(i,SAVE_CHK)+1:i,:), [svfld sprintf('etaStar%06d.fits',floor(i/SAVE_CHK)+1)]);
    fitswrite(iPSF_BB, [svfld sprintf('psfBB%06d.fits',floor(i/SAVE_CHK)+1)]);
    if isZernWFR
        % save rms of final wfr used after reconstruction with modulation
        fitswrite(rmsrecons(i-mod(i,SAVE_CHK)+1:i), [svfld sprintf('rmsrecons%06d.fits',floor(i/SAVE_CHK)+1)]);
    end
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

if isFigsVisEnd
    h = findobj('type', 'figure');
    set(h, 'Visible', 'on')
end