%{
    Script for estimating the light loss due to reduced area on KPIC VFN

    The KPIC VFN mode applies a circular aperture mask onto the vortex mask to mitigate
      scattering form the edges of the mask since the square vortex we have is smaller than
      the pupil at the pupil-stage ("coronagraph") plane. However, this circular mask clips the
      edges of the daytime pupil by a decent amount and the clips the edges of the nightime 
      pupil by a bit. This script determines the loss expected due to clipping at nighttime.

    The script was relevant for estimating the correction needed to determine "true" VFN null 
    on KPIC by accounting for the losses due to nulling effect vs. due to area reduction.
%}

clear; close all; 
addpath(genpath(fullfile('..','VFNlib')));

%% Input parameters 

%-- Provide regular parameters
% Define smapling info
N = 2^12; % Size of computational grid (NxN samples) 
% NOTE: use -20 below so that we have extra room for circumscribed apertures
apRad = N/2-20; % Aperture radius in samples 

%-- Define parameters for falco MFT propagator
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 2.45e-3;    %[m] pupil size 

%% Generate the coordinate system

%-- Coordinates in the focal plane
coordsPP = generateCoordinates(N);% Creates NxN arrays with coordinates 
% Helpful for plotting: get axes coords in pupil radii
xvalsPP = coordsPP.xvals/apRad;
yvalsPP = coordsPP.yvals/apRad;

%-- Key values for getting scaling right using falco MFT propagator
% "pixel" size in pupil
coordsPP.dx  = DPup/(2*apRad);   % meters/sample

%% Create array with pupil function

PUPIL = makeKeckPupil( 2*apRad, N );
%PUPIL = makeCircularPupil( apRad, N );

%-- Get norm for coupling fractions (simple sum since using MFT propagator)
totalPower0 = sum(abs(PUPIL(:)));

figure(1); 
imagesc(xvalsPP,yvalsPP,PUPIL); 
axis image;
axis([-1 1 -1 1]);
title('Pupil');
colorbar; 
colormap(parula(256));
drawnow;

%% Simulate smaller pupil masks within that aperture

%-- Normalize RHO so that it represents a circumscribed circlular aperture
    % Since makeKeckPupil goes edge-to-edge with 2*apRad, normalizing by
    % apRad results in a very small number of pixels outside of the RHO <=1
    % unit circle. This produces a 0.05% error but it's easy enough to fix
    % so let's just renorm so that there's no error!
% Find the outermost pixel within PUPIL
outermostpoint = max(coordsPP.RHO(logical(PUPIL)));
RHO = coordsPP.RHO/outermostpoint;

%-- Compute a smaller pupil size
r1 = 11.7 / 12.67;   % (11.7 is the known KPIC VFN diam. 12.94 is average of ellipse sides for DS diameter)
% elements within the smaller aperture
inR1 = RHO <= r1;
% calculate amount of light within r1
totalPowerR1 = sum(abs(PUPIL(inR1)));

fprintf('r1 is %0.4fx full diam and has %0.4f%% of the full area\n', r1, totalPowerR1/totalPower0*100)

%-- Display the apertures
figure(); 
imagesc(PUPIL + double(RHO < 1) +  double(inR1))
axis image;
title('Overlayed Apertures');