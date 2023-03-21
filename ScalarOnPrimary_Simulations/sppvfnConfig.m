%% Input Parameters
% HOME PATH FOR THIS FILE - \VFN-Simulations\ScalarOnPrimary_Simulations

%-- Provide regular parameters
% Define sampling info
inpar.N = 2^10; % Size of computational grid (NxN samples) 
inpar.apRad = inpar.N/2-4; % Aperture radius in samples
inpar.apDia0 = 2 * inpar.apRad;

%-- Define wavelength info
inpar.lambda0 = 2.2e-6; %central wavelength
% fracBW = 0.2; % \Delta\lambda/\lambda
inpar.fracBW = 0.1818; %\Delta\lambda/\lambda

%***************************************************************************%
inpar.numWavelengths = 51; % number of discrete wavelengths
%***************************************************************************%

inpar.lambdas = getWavelengthVec(inpar.lambda0, inpar.fracBW, inpar.numWavelengths);% array of wavelengths (meters)

inpar.keckD = 10.949;
% inpar.lam0OverD = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

%-- Define charge of the vortex mask at each wavelength
inpar.charge = 1;%*ones(1,inpar.numWavelengths); % achromatic

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
lambda0Fnum_samp = 51; %[samples/ (lam0 F#)] in units of samples in the image plane
im_size = 10;    %[lam0/D] field of view in final image plane

inpar.numRings = 3;
inpar.wGap = 25.4/10916*inpar.apDia0/2;