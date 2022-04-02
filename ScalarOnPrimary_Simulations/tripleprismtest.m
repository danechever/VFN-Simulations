close all;
clear;
clc;

%% Input parameters 

% Define smapling info
inPar.N = 2^11; % Size of computational grid (NxN samples) 
inPar.apRad = 128; % Aperture radius in samples 
inPar.apDia0 = 2 * inPar.apRad;

% Define wavelength info
inPar.lambda0 = 2.2e-6; %central wavelength
inPar.fracBW = 0.1818; %\Delta\lambda/\lambda
inPar.numWavelengths = 5;% number of discrete wavelengths 
inPar.lambdas = getWavelengthVec(inPar.lambda0,inPar.fracBW,inPar.numWavelengths);% array of wavelengths (meters)

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
inPar.charge = 1*(inPar.lambda0./inPar.lambdas); % simple scalar model

% Define wavefront error at the central wavelength
inPar.nolls = 4:8;
inPar.coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
inPar.offsetX = 0;%0.0952*apRad;
inPar.offsetY = 0;%0.0524*apRad; 

inPar.numRings = 3;
inPar.wGap = 25.4/10916*inPar.apDia0/2;

wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

inPar.wedge_ADC = 45 * pi/180; %Thorlabs 30 arcminute CaF2 wedge

n0 = 1;
n1 = getRefractiveIndex(wedge_ADC_1, 1e6*inPar.lambda0);
n2 = getRefractiveIndex(wedge_ADC_2, 1e6*inPar.lambda0);
n3 = getRefractiveIndex(wedge_ADC_3, 1e6*inPar.lambda0);

phi1 = 7.0516 * pi/180;
phi2 = 3.8050 * pi/180;
phi3 = 1.1465 * pi/180;

clocking = [1,1,1,1,1,1]*-90*pi/180;

tilt1 = 0;

dz = 0; %2*pi*(1.3751 * (inPar.lambdas/inPar.lambda0) - 1.3501);

    % angle going through first prism
    
u0x = sin(dz);
u0y = 0;
u0z = cos(dz);
u0 = [u0x, u0y, u0z];
    
norm1 = [0, 0, -1];
u1 = snell_in3d(u0, norm1, n0, n1);

% angle going into 2nd prism
phitot = (phi1);
norm2 = [(sin(phitot) * cos(clocking(1))), (sin(phitot) * sin(clocking(1))), -cos(phitot)];
u2 = snell_in3d(u1, norm2, n1, n2);
    
    % angle going into 3rd prism
phitot = (phi1 - phi2);    
norm3 = [(sin(phitot) * cos(clocking(2))), (sin(phitot) * sin(clocking(2))), -cos(phitot)];
u3 = snell_in3d(u2, norm3, n2, n3);
    
    % angle going into the air between triplets
phitot = (phi1 - phi2 -phi3);
norm3p = [(sin(phitot) * cos(clocking(3))), (sin(phitot) * sin(clocking(3))), -cos(phitot)];
u3p = snell_in3d(u3, norm3p, n3, n0);
    
    % angle going into 4th prism
phitot = (phi3 + phi2 - phi1);
norm4 = [(sin(phitot) * cos(clocking(4))), (-1*sin(phitot) * sin(clocking(4))), -cos(phitot)]; % clock other way now
u4 = snell_in3d(u3p, norm4, n0, n3);
    
    % angle going into 5th prism
phitot = (phi2 - phi1);
norm5 = [(sin(phitot) * cos(clocking(5))), (-1*sin(phitot) * sin(clocking(5))), -cos(phitot)]; % clock other way now
u5 = snell_in3d(u4, norm5, n3, n2);
    
    % angle going through 6th prism
phitot = (-phi1);
norm6 = [(sin(phitot) * cos(clocking(6))), (-1*sin(phitot) * sin(clocking(6))), -cos(phitot)]; % clock other way now
u6 = snell_in3d(u5, norm6, n2, n1);
    
    % angle leaving 6th prism relative to surface normal (which is also direction of chief ray)
phitot = 0;
norm6p = [0, 0, -1];
u6p = snell_in3d(u6, norm6p, n1, n0);
    
dz_OUT = atan2(u6p(1), u6p(3));



%dz_out = ADC_PRISM(dz, n0, n1_ADC, n2_ADC, n3_ADC, phi1_ADC, phi2_ADC, phi3_ADC, clocking);

