addpath(genpath(fullfile('..','VFNlib')));
addpath(genpath(fullfile('..','..','falco-matlab')))

%%
% First pass at making a script that will generate realistic CRED2 images
%    with pre-specified TT and WFE residuals

%-- First create the simulation structure
Npup = 2^7;     % pupil-plane sample grid size
Nim  = 64;      % focal-plane sample grid size
charge = 0;     % vortex charge to simulate
lam0 = 1640e-9; % central wavelength 
starHmag = 4.26;% H-mag of the star under consideration
% Create struct
cred2sim = setUp_CRED2ImSim(Npup, Nim, charge, lam0, starHmag);


%-- Run simulation to generate cube of data
nframes = 1000;
Tint = 0.028999; 
% Add semi-realistic WFE
    % Add a low-order trefoil with noll 9. Then add high-order stuff to get
    % the static speckle halo we get with the SHWFS on Keck
nolls = [9, 50:150];
    % Provide a trefoil coeff matched by eye. For the speckle halo, sample
    % from the standard uniform distribution with a peak amplitude set to
    % provide roughly the right counts/Strehl in the image (determined by
    % trial and error for now).
coeffs = 0*[0.08, 0.035*rand(1,length(nolls)-1)];  % [Waves RMS]
% Sample tilts from 2 gaussian distributions with different parameters
sig_x = 45;     % [mas] rms X-jitter 
sig_y = 45;     % [mas] rms Y-jitter
mu_x  = 0;      % [mas] average X position
mu_y  = 0;      % [mas] average Y position
tilts = [normrnd(mu_x, sig_x, [nframes,1]), normrnd(mu_y,sig_y,[nframes,1])];
% Decide whether iteration counter should be printed
verbose = true;
% Generate images
[noisy, ~] = CRED2ImSim_getNoisyIm(cred2sim, Tint, nolls, coeffs, tilts, verbose);

%-- OPTIONALS:
% Option to display frames. NOTE: if more than 30 frames are generated,
    % then only 30 are displayed, evenly sampled in data cube.
isDispFrames = true;
% Option to save fits file at end as well as filename to use
isSaveFits = false;
pth = '/media/Data_Drive/KPIC/dev/dechever/TTJitCharac/Algo_Dev/CRED2_Sim_SampleData/';
flnm = sprintf('%05dSamp_%06.3fWFECoeff_TiltMean%05.2f%05.2f_RMS%05.2f%05.2f.fits',nframes,frstcoeff,mu_x,mu_y,sig_x,sig_y);
flnm = [pth flnm];

%% OPTIONAL: Display the images
if isDispFrames
    % Make sure we only display up to 30 frames
    NFrames = max([size(coeffs,1),size(tilts,1)]);
    if NFrames > 30
        % If more than 30 frames would be displayed, evenly sample 30
        step = ceil(NFrames/30);
    else
        step = 1;
    end
    % Define counter such that figure handle numbers are still sequential
    fig_ctr = 1;
    for frame = 1:step:NFrames
        figure(fig_ctr+200);
        fig_ctr = fig_ctr + 1;
        imagesc(squeeze(noisy(frame, :,:)));
        % add crosshair on commanded tilt position
        xt = tilts(frame,1)/(cred2sim.lambda0OverKeckD_mas/cred2sim.lambda0Fnum_samp) + Nim/2 + 1;
        yt = tilts(frame,2)/(cred2sim.lambda0OverKeckD_mas/cred2sim.lambda0Fnum_samp) + Nim/2 + 1;
        hold on; scatter(xt, yt, 40, 'r+'); hold off;
        title(sprintf('Frame #%d',frame));
        xlabel('[pix]');
        ylabel('[pix]');
        %colormap('gray');
        colorbar;
        drawnow
    end
end

%% OPTIONAl: Save resulting cube into a fits file
if isSaveFits
    % First check if the provided filename is unique.
    if isfile(flnm)
        % file exists so warn user
        error('Filename already exists. You can provide new name in "flnm" variable and then re-run just the fits saving section of the script')
    elseif ~isfolder(fileparts(flnm))
        % filename is unique but requested directory doesn't exist so make
        % it for the user
        warning('Requested directory for file did not exist, we will make it for you')
        mkdir(fileparts(flnm));
    end
    import matlab.io.*
    fitmap= fits.createFile(flnm);
    % Add the noisy images first
        % NOTE: permute is needed in order to get X,Y correct in python import
    fits.createImg(fitmap,'double',size(permute(noisy,[3,1,2])));
    %header data
    fits.writeKey(fitmap, 'NAX1', 'frames');
    fits.writeKey(fitmap, 'NAX2', 'X-axis');
    fits.writeKey(fitmap, 'NAX3', 'Y-axis');
    fits.writeKey(fitmap, 'Npup', Npup);
    fits.writeKey(fitmap, 'Nim', Nim);
    fits.writeKey(fitmap, 'charge', charge);
    fits.writeKey(fitmap, 'lam0', lam0);
    fits.writeKey(fitmap, 'starHmag', starHmag);
    fits.writeKey(fitmap, 'tint', Tint);
    fits.writeKey(fitmap, 'sigX', sig_x);
    fits.writeKey(fitmap, 'sigY', sig_y);
    fits.writeKey(fitmap, 'muX', mu_x);
    fits.writeKey(fitmap, 'muY', mu_y);
    fits.writeKey(fitmap, 'apRad', cred2sim.apRad);
    fits.writeKey(fitmap, 'DPup', cred2sim.DPup);
    fits.writeKey(fitmap, 'foc', cred2sim.foc);
    fits.writeKey(fitmap, 'KeckD', cred2sim.Keck_D);
    fits.writeKey(fitmap, 'lam0Dmas', cred2sim.lambda0OverKeckD_mas);
    fits.writeKey(fitmap, 'lam0Dpix', cred2sim.lambda0Fnum_samp);
    fits.writeKey(fitmap, 'CRED2QE', cred2sim.CRED2_QE);
    fits.writeKey(fitmap, 'CRED2DC', cred2sim.CRED2_DarkCurrent);
    fits.writeKey(fitmap, 'CRED2RN', cred2sim.CRED2_ReadNoiseRMS);
    fits.writeKey(fitmap, 'CRED2ADU', cred2sim.CRED2_ADU);
    fits.writeKey(fitmap, 'E2ETHPT', cred2sim.KPIC_thpt);
    fits.setCompressionType(fitmap,'NOCOMPRESS');
    fits.writeImg(fitmap,permute(noisy,[3,1,2]));
    % Add the TT and WFE values
    fits.createImg(fitmap, 'double', size(tilts));
    fits.writeImg(fitmap,tilts);
    fits.createImg(fitmap, 'int16', size(nolls));
    fits.writeImg(fitmap,nolls);
    fits.createImg(fitmap, 'double', size(coeffs));
    fits.writeImg(fitmap,coeffs);
    fits.closeFile(fitmap);
end