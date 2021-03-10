function mask = generateVortexMaskKeckPrimary(inputs, seed)
%Function to generate the vortex mask applied to the keck pupil through
%applying tip/tilt/piston to the segmented mirrors
%% Input parameters 

% Define smapling info
param.N = inputs.N; % Size of computational grid (NxN samples) 
param.apRad = inputs.apRad; % Aperture radius in samples 
param.apDia0 = inputs.apDia0;

% Define wavelength info
param.lambda0 = inputs.lambda0; %central wavelength
param.fracBW = inputs.fracBW; %\Delta\lambda/\lambda
param.numWavelengths = inputs.numWavelengths;% number of discrete wavelengths 
param.lambdas = inputs.lambdas;% array of wavelengths (meters)

%Define charge of the vortex mask at each wavelength
%charge = ones(1,numWavelengths); % achromatic
param.charge = inputs.charge; % simple scalar model

% Define wavefront error at the central wavelength
param.nolls = inputs.nolls;
param.coeffs = inputs.coeffs;

% Give offsets for the vortex mask
param.ofssetX = inputs.offsetX;%0.0952*apRad;
param.offsetY = inputs.offsetY;%0.0524*apRad; 

param.numRings = inputs.numRings;
param.wGap = inputs.wGap;

%% Generate the coordinate system

%coords = generateCoordinates(Sinf.N);% Creates NxN arrays with coordinates 
xvals = inputs.xvals;% Helpful for plotting
yvals = inputs.yvals;

%% Create matrix with tilt and pistons for the analytical solution, and
%%generate constant coefficients for phase and amplitude for each hexagonal segment.

hexFlatDiam = (param.apDia0-3*2*param.wGap)/(2*3+1);
param.hexSep = hexFlatDiam + param.wGap;

if (~exist('seed', 'var'))
   param.hexAmpConst = NaN(param.N, param.N, 36);
   param.hexPhzConst = NaN(param.N, param.N, 36);
   
   initial = NaN(36,3); %xtilt, ytilt, piston

   loc = 1;
   for ringNum = 1:3 %adds absolute value of tilt to each segment in the position they are added.
       for seg = (loc):(loc+ringNum*6-1)
           initial(seg,1:2) = (1*2*pi)/(ringNum*6);
       end
       loc = (loc+ringNum*6);
   end
   
   segs = ones(1,36);
   count = 1;
   for ringNum = 1:3
       crow = ringNum * param.hexSep;
       ccol = 0;
       t = atan2(crow,ccol);
       t = round(t,3);
       
       if(segs(count) == 1)
           %Uncomment below to use analytical solution / Comment out the 3
           %lines below if using a seed
           initial(count,1) = initial(count,1) * -sin(t);
           initial(count,2) = initial(count,2) * cos(t);
           initial(count,3) = t/(2*pi);
           
           [param.hexAmpConst(:,:,count), param.hexPhzConst(:,:,count)] = ...
               generateHexConstants(crow, ccol, param.numRings, param.apDia0, param.wGap, zeros(param.N));
       end
       
       count = count + 1;
       for face = 1:6
           step_dir = pi/6*(2*face+5);
           steprow = param.hexSep*sin(step_dir);
           stepcol = param.hexSep*cos(step_dir);
           stepnum = 1;
           
           while(stepnum <= ringNum && count <= 36)
               crow = crow + steprow;
               ccol = ccol + stepcol;
               
               t = atan2(crow,ccol);
               t = round(t,3);
               
               if(face==6 && stepnum ==ringNum)
                   disp('finished ring');
               else
                   if(segs(count) == 1)
                       %Uncomment below to use analytical solution / Comment
                       %out the 3 lines below if using a seed
                       initial(count,1) = initial(count,1) * -sin(t);
                       initial(count,2) = initial(count,2) * cos(t);
                       initial(count,3) = t/(2*pi);
                       
                       %consider third output that is all 5 different wedge
                       %effects
                       [param.hexAmpConst(:,:,count), param.hexPhzConst(:,:,count)] = ...
                           generateHexConstants(crow, ccol, param.numRings, param.apDia0, param.wGap, zeros(param.N));
                   end
                   count = count + 1;
               end
               stepnum = stepnum + 1;
           end
       end
   end
else
    initial = seed;
end
 
%% Create array with pupil function
addpath(['..' filesep '..' filesep 'falco-matlab' filesep 'lib' filesep 'utils']);
%-- Decrease matrix size in pupil plane to reduce runtime
param.PUP_CRP_SZ = round(2.1*param.apRad);
param.hexAmpConst = pad_crop(param.hexAmpConst,param.PUP_CRP_SZ);
param.hexPhzConst = pad_crop(param.hexPhzConst,param.PUP_CRP_SZ);

%% Define pupil field

mask = angle(makeKeckPupilInputs(param, initial));
