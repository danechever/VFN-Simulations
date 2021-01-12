%{
Code for summarizing key results of an on-sky null simulation

*** This script analyzes data from KPICVFN_OnSky_NullCalculator
*** Output of this script can be used as input for RMS2Null_Plotter

Goals:___
- Take in null depth and RMS of the reconstructed WF from a simulation
- Report the Mean and STD of those parameters
- Determine rough approximation bsaed on inversion and user input
  * NOTE: the RMS2Null_Plotter script is better for getting the approx.

Notes:___
* This code only analyzes a single simulation (ie. a single output of the
    NullCalculator code so any approx. from this will be based on a single
    point. This resulted in a good approximation in general but the
    RMS2Null_Plotter script does a more in-depth and robust approx. which
    should be more reliable. 
%}

clear; close all; 
addpath(genpath(['..' filesep 'VFNlib']));
addpath(genpath(['..' filesep '..' filesep 'VFN-Lab' filesep 'AnalysisCode']))
addpath(genpath(['..' filesep '..' filesep 'falco-matlab']))

%% Input parameters 

%-- Extract values below from fits header
% Folder with files
wfrfld = '/media/Data_Drive/VFN/KPIC_PyWFS_Data/telemetry_20190617/Processed/10msSampling/';
datfld = 'ZernAllx001p000_Charge1/';
% SET exponent for approximations
pwr = 2;
% SET COEFF for validation
coeff = 1.51;
% filename counter format
wfrFMT = '%06d';
% fits to reference
fitsnm = sprintf(['wfr_reform' wfrFMT '.fits'],1);
% Extract keywords
kwds    = fitsinfo([wfrfld fitsnm]);
kwds    = kwds.PrimaryData.Keywords;      % all keywords

%-- Get total number of samples from keywords
wfSamps = VFN_An_getKwd(kwds, 'wfsamps');  % Old method: size(dir([wfrfld wfrfnm '*.fits']),1);

%% Read in data

%-- Preallocate rms  and null matrices
% Get keywords which describe shape of files
nullnm = 'etaStar';
rmsfnm = 'rmsrecons';
kwdsN    = fitsinfo([wfrfld datfld sprintf([nullnm wfrFMT '.fits'],1)]);
kwdsN    = kwdsN.PrimaryData.Keywords;      % all keywords
frames_per_file = VFN_An_getKwd(kwdsN, 'NAXIS2');
% Preallocate full file
numWavelengths = VFN_An_getKwd(kwdsN, 'NAXIS1');
nulls = nan(wfSamps, numWavelengths); 
rmsWF = nan(wfSamps, 1);        % use 1 since we know rms matrix has 1 col

%-- Read all files and place in null/rmsWF matrix
% Get total number of files over which to iterate
N_fls = size(dir([wfrfld datfld nullnm '*.fits']),1);
for i = 1:N_fls
    % read the files
    null_tmp = fitsread([wfrfld datfld sprintf([nullnm wfrFMT '.fits'],i)]);
    rms_tmp = fitsread([wfrfld datfld sprintf([rmsfnm wfrFMT '.fits'],i)]);
    % Place in matrix
    if i == N_fls
        % On the last file so change indexing since it may have less frames
        nulls((frames_per_file*(i-1))+1:end, :) = null_tmp;
        rmsWF((frames_per_file*(i-1))+1:end, :) = rms_tmp;
    else
        % On every other file, simply place it in the matrix
        nulls((frames_per_file*(i-1))+1:frames_per_file*i, :) = null_tmp;
        rmsWF((frames_per_file*(i-1))+1:frames_per_file*i, :) = rms_tmp;
    end
end  

%-- Check if any frames were left unallocated
if any(isnan(nulls),'all')
    error('A zernike coeff value was left unallocated')
end

clear rms_tmp N_fls i null_tmp

%% Calculate coefficient
cent_lam = ceil(numWavelengths/2);

fprintf('----> SAMPLE: %s\n', datfld);
%-- Single (central) wavelength
%- Calculate means then coeff
mn_null = mean(nulls(:,cent_lam));
mn_rms = mean(rmsWF);
% Calculate coeff
cent_coeff_mn1st = (mn_null^(1/(pwr))) / mn_rms;
fprintf('Cent lam mean first: %f\n', cent_coeff_mn1st)

%- Calculate coeff then mean
%calculate coeff
cent_coeff_mn2nd = (nulls(:,cent_lam).^(1/(pwr))) ./ rmsWF;
% Calculate mean
cent_coeff_mn2nd = mean(cent_coeff_mn2nd);
fprintf('Cent lam mean second: %f\n', cent_coeff_mn2nd)

%- Validate
fprintf('Ideal Approx Eq: null = (%0.3f*rms)^(%d)\n', cent_coeff_mn1st, pwr)
% Calculate null using approximation
approx = (cent_coeff_mn1st*rmsWF).^(pwr);
bnd = 5;    % value within which approx should be from null
val = sum((approx < bnd*nulls(:,cent_lam)) & (approx > nulls(:,cent_lam)/bnd));
fprintf('Num times approx is within a factor of %0.1f of null:\n', bnd);
fprintf('    %d out of %d\n', val, length(approx));

%- Validate with provided number
fprintf('Manual Approx Eq: null = (%0.3f*rms)^(%d)\n', coeff, pwr)
% Calculate null using approximation
approx2 = (coeff*rmsWF).^(pwr);
bnd = 5;    % value within which approx should be from null
val2 = sum((approx2 < bnd*nulls(:,cent_lam)) & (approx2 > nulls(:,cent_lam)/bnd));
fprintf('Num times man approx is within a factor of %0.1f of null:\n', bnd);
fprintf('    %d out of %d\n', val2, length(approx2));

%-- Old error reporting with issues
% % Calculate fractional error on approximation
% err = (approx-nulls(:,cent_lam))./nulls(:,cent_lam);
% % display median error
% fprintf('Median fractional error: %f\n', median(err));
% display number of instances where approx is within factor of 2 of true val
% OLD, wrong way of doing bound: sum(abs(err)<2); 

%-- Display rmsWFE and null statistics
% compute STDs
std_rms = std(rmsWF);
std_null = std(nulls(:,cent_lam));
fprintf('- Mean WFE:   %f\n', mn_rms)
fprintf('- STD WFE:    %f\n', std_rms)
fprintf('- Mean Null:  %0.8f\n', mn_null)
fprintf('- STD Null:   %0.8f\n', std_null)

fprintf('\n\n')

