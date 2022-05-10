% Demo script showing how to use the CRED2ImSim functions
close all;

%-- First create the simulation structure
Npup = 2^7;     % pupil-plane sample grid size
Nim  = 64;      % focal-plane sample grid size
charge = 0;     % vortex charge to simulate
lam0 = 1640e-9; % central wavelength 
starHmag = 5.26;% H-mag of the star under consideration
% Create struct
cred2sim = setUp_CRED2ImSim(Npup, Nim, charge, lam0, starHmag);


%-- Run simulation for 9 frames with no WFE but unique TT error
Tint = 0.003;
nolls = [5,8,9];
coeffs = 0*[-0.08, 0.05, 0.09];
nframes = 9;
% Sample tilts from 2 gaussian distributions with different parameters
sig_x = 10;     % [mas] rms X-jitter 
sig_y = 10;     % [mas] rms Y-jitter
mu_x  = 0;      % [mas] average X position
mu_y  = 0;      % [mas] average Y position
tilts = [normrnd(mu_x, sig_x, [nframes,1]), normrnd(mu_y,sig_y,[nframes,1])];
% Generate images
[noisy, psf] = CRED2ImSim_getNoisyIm(cred2sim, Tint, nolls, coeffs, tilts);

%-- Display results
figure();
for im = 1:size(noisy,1)
    subplot(3,3,im);
    imagesc(squeeze(noisy(im,:,:)))
    axis image
end