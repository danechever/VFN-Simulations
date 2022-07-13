%% ADC Tilt Optimizer Function
% This function is designed to determine the optimal clocking angle of the 
% ADC for correcting VFN introduced null-shift.
% Working correctly as of last edit 

clear; close all;
addpath(['..' filesep 'VFNlib']);
load pyyo;
load oldetas;

% py.importlib.import_module('sellmeir1');
% py.importlib.import_module('triple_prism');
sellmeir1coeffs;
%% Input parameters 

% Define smapling info
inpar.N = 2^11; % Size of computational grid (NxN samples) 
inpar.apRad = 128; % Aperture radius in samples 
inpar.apDia0 = 2 * inpar.apRad;
inpar.keckD = 10.949; %Meters

% Define wavelength info
inpar.lambda0 = 2.2e-6; %central wavelength
inpar.fracBW = 0.1818; %\Delta\lambda/\lambda
inpar.numWavelengths = 5;% number of discrete wavelengths 
inpar.lambdas = getWavelengthVec(inpar.lambda0,inpar.fracBW,inpar.numWavelengths);% array of wavelengths (meters)
inpar.lam0OverD = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
inpar.charge = 1*(inpar.lambda0./inpar.lambdas); % simple scalar model

% Define wavefront error at the central wavelength
inpar.nolls = 4:8;
inpar.coeffs = 0*[0.01,-0.0099,-0.0095,-0.0008,0.0033];

% Give offsets for the vortex mask
inpar.offsetX = 0;%0.0952*apRad;
inpar.offsetY = 0;%0.0524*apRad; 

inpar.numRings = 3;
inpar.wGap = 25.4/10916*inpar.apDia0/2;

lin_offset = [1.1001 -1.1001];
inpar.p = lin_offset(1);

wedge_ADC_1 = 'BaF2';
wedge_ADC_2 = 'CaF2';
wedge_ADC_3 = 'ZnSe';

%Wedge angle of Single Wedge Prism correcting wavelength dependent offset
sprism_angle = atan(890.16*(inpar.p*inpar.lambda0)/(inpar.keckD*(getRefractiveIndex(wedge_ADC_2 ,1e6*inpar.lambda0) - 1)));

tilt_vals = 0;

n0 = zeros(1,inpar.numWavelengths) + 1;

%ADC Wedge Angles
phi1 = 7.0516 * pi/180; 
phi2 = 3.8050 * pi/180;
phi3 = 1.1465 * pi/180;

%offset_lamoverd = [-0.1000   -0.0375    0.0250    0.0875    0.1500];

%Clocking angle (2x this is the actual angle between triplets)
clocking = -34.1697;%-89.988;

%dz = zeros(1,5);

%Loop computing input deviation in radians and output deviation for a given
%clocking angle

wvs = py.numpy.array(1e6*inpar.lambdas);
n1 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));
n2 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));
n3 = py.sellmeir1.sellmeir1(wvs, 273.15 + 3, znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));

for ch = 1:inpar.numWavelengths
    
    %disp(lin_offset(1)*(inpar.lambdas(ch)/inpar.lambda0) + lin_offset(2));
    
    %dz(ch) = inpar.lambda0/inpar.keckD*(lin_offset(1)*(inpar.lambdas(ch)/inpar.lambda0) + lin_offset(2));
    %dz = 0;
    %disp(dz(ch));
    
    dz(ch) = (-1/890.16)*(getWedgeTilt(wedge_ADC_2, sprism_angle, 1e6*inpar.lambdas(ch))...
        - getWedgeTilt(wedge_ADC_2, sprism_angle, 1e6*inpar.lambdas(ceil(inpar.numWavelengths/2))));

%     n1(ch) = sellmeir1(2,0,baf2_args(1,:),baf2_args(2,:),baf2_args(3,:),baf2_args(4,:),baf2_args(5,:),baf2_args(6,:));%getRefractiveIndex(wedge_ADC_1, 1e6*inpar.lambdas(ch));
%     n2(ch) = sellmeir1(2,0,caf2_args(1,:),caf2_args(2,:),caf2_args(3,:),caf2_args(4,:),caf2_args(5,:),caf2_args(6,:));%getRefractiveIndex(wedge_ADC_2, 1e6*inpar.lambdas(ch));
%     n3(ch) = sellmeir1(2,0,znse_args(1,:),znse_args(2,:),znse_args(3,:),znse_args(4,:),znse_args(5,:),znse_args(6,:));%getRefractiveIndex(wedge_ADC_3, 1e6*inpar.lambdas(ch));
       
%     getADCTilt(dz(ch),1,n1(ch),n2(ch),n3(ch), phi1, phi2, phi3, clocking);
  
end

% dz = 1e-7*[-1.2119   -0.6059         0    0.6059    1.2119]; ???
% dz = 1e-7*[0   0        0    0    0];
dz = pyyCopy * inpar.lam0OverD;
disp(dz)
%% Show input dispersion and output dispersion for a given clocking angle (defined above)
disp(dz);

figure()
hold on
title(['Tilt across DS K Band from VFN on-sky']);
xlabel(['Wavelength (um)']);
ylabel(['Tilt (mas)']);
plot(1e6* inpar.lambdas,dz * (3600*1000 *180)/pi,'Color', 'b');
txt = ['dispersion across band: ' num2str((3600*1000 *180)/pi*(dz(end) - dz(1)))];
text(2.2,0,txt);
hold off

% dz = [0 0 0 0 0 0 0 0 0 0 0] ;
dz = 890.16*dz;
dz = py.tuple(dz);
dz = py.numpy.array(dz);
tilt = py.triple_prism.triple_prism(dz,n0,n1,n2,n3,phi1,phi2,phi3,clocking,tilt_vals);
    
tilt = tilt.tolist();

for j = 1:inpar.numWavelengths
    tilt_out(j) = tilt{j};
end

tilt_out = tilt_out / 890.16;

tilt_out = tilt_out - tilt_out(ceil(inpar.numWavelengths / 2));

figure()
hold on
title(['Tilt across DS K Band @ 34.1697 deg on-sky']);
xlabel(['Wavelength (um)']);
ylabel(['Tilt (mas)']);
plot(1e6* inpar.lambdas,tilt_out * (3600*1000 *180)/pi,'Color', 'b');
txt = ['dispersion across band: ' num2str((3600*1000 *180)/pi*(tilt_out(end) - tilt_out(1)))];
text(2.2,0,0,txt);
hold off

dz1 = zeros(1,inpar.numWavelengths); %[0 0 0 0 0 0 0 0 0 0 0] ;
dz1 = py.tuple(dz1);
dz1 = py.numpy.array(dz1);
tilt_vals = py.triple_prism.triple_prism(dz1,n0,n1,n2,n3,phi1,phi2,phi3,clocking,tilt_vals);
    
tilt_vals = tilt_vals.tolist();

for j = 1:inpar.numWavelengths
    tilt_out1(j) = tilt_vals{j};
end

tilt_out1 = tilt_out1 - tilt_out1(ceil(inpar.numWavelengths / 2));

% n1 = py.numpy.array(n1);
% n2 = py.numpy.array(n2);
% n3 = py.numpy.array(n3);
% n0 = py.numpy.array(n0);
% dz = py.numpy.array(dz);


%tilt = py.triple_prism.triple_prism(dz,n0,n1,n2,n3,phi1,phi2,phi3,clocking,tilt1);


% global_offset = 1*getADCTilt(dz(3),1,n1(3),n2(3),n3(3), phi1, phi2, phi3, clocking);

%tilt1 = getADCTiltOpt(dz,n0,n1,n2,n3,phi1,phi2,phi3,clocking);

% disp(dz);
% disp(tilt - global_offset);
%disp(tilt1);

%disp(dz*206265*1000);
%disp(tilt*206265*1000);

%% Gridscan optimization

%-- Parameters for SMF (Thorlabs SM600 in this case)
%     % link: https://www.thorlabs.com/NewGroupPage9_PF.cfm?ObjectGroup_ID=949
fiber_props.core_rad = 5.5e-6;% Core radius [um]
fiber_props.n_core = 1.4436;% core index (interp @220nm)
fiber_props.n_clad = 1.4381;% cladding index (interp @220nm)
fiber_props.type = 'gaussian';

%-- Define parameters for falco MFT propagator
%<> NEED TO FIGURE THESE OUT
    % these don't matter in themselves as long as they are consistent w/ each 
    % other and with lambda
% OPTION 1: define focal length and pup diameter manually
%foc = 11e-3;      %[m] final focal length
DPup = 2.45e-3;    %[m] pupil size 
% OPTION 2: solve for focal length based on ideal Fnum and pup diameter
Fnum = getMFD(fiber_props,inpar.lambda0)/(inpar.lambda0*1.4); % focal ratio of the beam at the fiber
foc = Fnum*DPup;

lambda0Fnum_meters = inpar.lambda0*foc/DPup;

inpar.lam0OverD_meters = lambda0Fnum_meters/foc; %MFT value

%%
dz = pyyCopy*inpar.lam0OverD;

% dz = pyyCopy*inpar.lam0OverD_meters; %MFT Value

dz = dz*890.16;
dz = py.tuple(dz);
dz = py.numpy.array(dz);

samp = linspace(0,90,100);
tilt_vals = [];

for i = 1:100
    
    clocking = deg2rad(samp(i));
%     disp(clocking);
    tilt = py.triple_prism.triple_prism(dz,n0,n1,n2,n3,phi1,phi2,phi3,clocking,0);
    
    tilt = tilt.tolist();
    
    for j = 1:inpar.numWavelengths
        tilt_out(j) = tilt{j};
    end
    
    tilt_out = tilt_out - tilt_out(ceil(inpar.numWavelengths/2));
    %disp(tilt_out);
    tilt_out = tilt_out/890.16;
    
    % Wavelength dependent dispersion
    tilt_vals(:,i) = tilt_out';
%     for j = 1:11
%         tilt_vals(:,i) = tilt_out;
%     end
    
    % Sum of dispersion at each wavelength squared
    tiltsums(i) = sum(tilt_out.^2);
    
    %Slope of ADC induced tilt
    tiltslopes(i) = tilt_out(end) - tilt_out(1);
    
end

slopes = tilt_vals(2,:) - tilt_vals(1,:);

figure()
hold on
title(['Gridsearch Min']);
xlabel(['Clocking Angle 0-90 deg']);
ylabel(['Sum of Squares']);
% plot(samp, tilt1,'Color', 'r');
% plot(samp, tilt2,'Color', 'y');
% plot(samp, tilt3,'Color', 'g');
% plot(samp, tilt4,'Color', 'c');
% plot(samp, tilt5,'Color', 'm');
plot(samp, tiltsums, 'Color', 'b');
hold off

[M,I] = min(tiltsums);
disp(M);
disp(I);
disp(samp(I));
disp(tilt_vals(:,I));

figure()
hold on
title(['Polychromatic Tilt 2 - 2.4 um']);
xlabel(['Clocking Angle 0 - 90 deg']);
ylabel(['Tilt (radians)']);
plot(samp, tilt_vals,'Color','r');

for i = 1:inpar.numWavelengths
    plot(samp, tilt_vals(i,:));
end

% plot(samp, tilt_vals,'Color', 'r');
% plot(samp, tilt2,'Color', 'y');
% plot(samp, tilt3,'Color', 'g');
% plot(samp, tilt4,'Color', 'c');
% plot(samp, tilt5,'Color', 'm');
% plot(samp, tilt6,'Color', 'r');
% plot(samp, tilt7,'Color', 'y');
% plot(samp, tilt8,'Color', 'g');
% plot(samp, tilt9,'Color', 'c');
% plot(samp, tilt10,'Color', 'm');
% plot(samp, tilt11,'Color', 'r');
% plot(samp, tiltsums, 'Color', 'b');
hold off

figure()
hold on
title(['Slope of ADC tilt vs increasing clocking angle']);
xlabel(['Clocking Angle 0 - 90 deg']);
ylabel(['ADC Dispersion Slope (tilt @ 2.1um - tilt @ 2.0um)']);
plot(samp, slopes,'Color','b');
hold off

figure()
hold on
title(['Dispersion across band for clocking angles 0-90deg']);
xlabel(['Clocking Angle 0 - 90 deg']);
ylabel(['Dispersion across band (mas) (tilt @ 2.4um - tilt @ 2.0um)']);

plot(samp, tiltslopes * (3600*1000 *180)/pi,'Color','b');

hold off
%% Deviation from normal of ADC output surface minimizing function
%nItrCap = 1000;
%vars.relTintArray = nan(nItrCap*1.2,1);
%vars.n_calls = 1;
% devdifMin = @(clocking) optClock(dz, n0, n1, n2, n3, phi1, phi2, phi3, clocking,tilt1);

%% Optimization
%options = optimoptions('fmincon','MaxFunctionEvaluations');
% [localMinX,fvalY] = fminsearch(devdifMin, clocking);
% disp(localMinX);
% disp(fvalY);






