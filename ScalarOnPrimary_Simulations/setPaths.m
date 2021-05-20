% setPaths.m

% Add paths to FALCO and PROPER 
mp.path.falco = '/home/vfndev/Documents/MATLAB/falco-matlab'; 
mp.path.proper = '/home/vfndev/Documents/MATLAB/SupportPackages/proper_v3.2_matlab_11feb20';

addpath(genpath(mp.path.falco),mp.path.proper);
