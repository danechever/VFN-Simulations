%{
Code for determining RMS2Null approximation via Logarithmic Fit

*** This script analyzes data from KPICVFN_OnSky_NullCalculator
*** RM2Null_CoeffCalculator outputs the values needed here

Goals:___
- Take in mean and STD of null depth and RMS WF from simulations
- Determine the ideal approximation/model to map WFE to Null depth
- Allow for manual tuning of the approximation (ex: check quadratic form
    explicitly vs. just using whatever power provides the ideal fit)
- Plot the results and allow for easy visualization
- Plot the zernike coefficients used in the simulation (almost like a PSD)

Notes:___
* As written here, the code has manual inputs for the mean and std values
    needed in the approximation. These values were taken from the results
    printed by RMS2Null_CoeffCalculator. This was done for simplicity on
    the day I made this code.
    - However, it should be easy to make the CoeffCalc. code save the
      results somewhere and read them in here if someone wants to put the
      time in to make those mods :)
%}


close all; clear all;

%% RMS to Null Analysis

%-- Populate data matrix (using data from VFN server KPIC simulations)
%########### CHARGE 2 ####################
% Layout of matrix (row-wise)
% NullMult, WFE mean, WFE STD, Null mean, Null STD
dat_ch2 = [100, 4.113753, 1.141707, 0.000278,   0.000364;
        50, 2.056876, 0.570854, 0.001086,   0.001332;
        25, 1.028438, 0.285427, 0.004156,   0.004893;
        10, 0.411375, 0.114171, 0.020017,   0.019226;
         5, 0.205688, 0.057085, 0.035043,   0.030588;
       2.5, 0.102844, 0.028543, 0.020624,   0.019833;
       1.0, 0.041138, 0.011417, 0.004494,   0.004862;
       0.5, 0.020569, 0.005709, 0.001183,   0.001371;
      0.25, 0.010284, 0.002854, 0.000299,   0.000355;
      0.10, 0.004114, 0.001142, 0.000048,   0.000057;
      0.05, 0.002057, 0.000571, 0.000012,   0.000014;
     0.025, 0.001028, 0.000285, 0.00000307, 0.00000361;
     0.010, 0.000411, 0.000114, 0.00000067, 0.00000071;
     0.005, 0.000206, 0.000057, 0.00000035, 0.00000028;
     0.003, 0.000123, 0.000034, 0.00000029, 0.00000016;
     0.001, 0.000041, 0.000011, 0.00000027, 0.00000005];
 
%-- Populate data matrix (using data from VFN server KPIC simulations)
%########### CHARGE 1 ####################
% Layout of matrix (row-wise)
% NullMult, WFE mean, WFE STD, Null mean, Null STD
dat_ch1 = [100, 4.113753, 1.141707, 0.00028249,   0.00037375;
        50, 2.056876, 0.570854, 0.00108922, 0.00140677;
        25, 1.028438, 0.285427, 0.00439101, 0.00535871;
        10, 0.411375, 0.114171, 0.02208851, 0.02127949;
         5, 0.205688, 0.057085, 0.02320178, 0.02240673;
       2.5, 0.102844, 0.028543, 0.00797345, 0.00940517;
       1.0, 0.041138, 0.011417, 0.00123715, 0.00156488;
       0.5, 0.020569, 0.005709, 0.00030083, 0.00036599;
      0.25, 0.010284, 0.002854, 0.00007475, 0.00008780;
      0.10, 0.004114, 0.001142, 0.00001226, 0.00001401;
      0.05, 0.002057, 0.000571, 0.00000332, 0.00000368;
     0.025, 0.001028, 0.000285, 0.00000105, 0.00000107;
     0.010, 0.000411, 0.000114, 0.00000039, 0.00000029;
     0.005, 0.000206, 0.000057, 0.00000029, 0.00000013;
     0.003, 0.000123, 0.000034, 0.00000027, 0.00000007;
     0.001, 0.000041, 0.000011, 0.00000025, 0.00000002];
 
%-- Flip the data so values increase downward
dat_ch1 = flipud(dat_ch1);
dat_ch2 = flipud(dat_ch2);

%-- Extract pertinent columns
% Charge 1
wfe_mn_ch1 = dat_ch1(:,2);
wfe_std_ch1 = dat_ch1(:,3);
nll_mn_ch1 = dat_ch1(:,4);
nll_std_ch1 = dat_ch1(:,5);
% Charge 2
wfe_mn_ch2 = dat_ch2(:,2);
wfe_std_ch2 = dat_ch2(:,3);
nll_mn_ch2 = dat_ch2(:,4);
nll_std_ch2 = dat_ch2(:,5);
 
%-- Perform log fit to data 
% Charge 1
ind2Fit_st_ch1 = 5;     %Index of point to start fit on
ind2Fit_nd_ch1 = 12;    %Index of point to end fit on
[P_ch1,S_ch1] = polyfit(log10(wfe_mn_ch1(ind2Fit_st_ch1:ind2Fit_nd_ch1)),log10(nll_mn_ch1(ind2Fit_st_ch1:ind2Fit_nd_ch1)),1);
xFit_ch1 = logspace(log10(min(wfe_mn_ch1)),log10(max(wfe_mn_ch1)),5);
yFit_ch1 = 10^P_ch1(end)*xFit_ch1.^P_ch1(1);
% Charge 2
ind2Fit_st_ch2 = 4;     %Index of point to start fit on
ind2Fit_nd_ch2 = 10;    %Index of point to end fit on
[P_ch2,S_ch2] = polyfit(log10(wfe_mn_ch2(ind2Fit_st_ch2:ind2Fit_nd_ch2)),log10(nll_mn_ch2(ind2Fit_st_ch2:ind2Fit_nd_ch2)),1);
xFit_ch2 = logspace(log10(min(wfe_mn_ch2)),log10(max(wfe_mn_ch2)),5);
yFit_ch2 = 10^P_ch2(end)*xFit_ch2.^P_ch2(1);


% Flag to show the manual approximations.
isFitANDApprox = true;
isApproxOnly = false;
%-- Set manual Approximations (Push towards quadratic fit)
manSlp_ch1 = 0.840;
yFitMan_ch1 = (manSlp_ch1*xFit_ch1).^2;
manSlp_ch2 = 1.65;
yFitMan_ch2 = (manSlp_ch2*xFit_ch2).^2;


% %-- SIMPLE Plot (no error bars)
% figure();
% loglog(wfe_mn,nll_mn, 'o', 'MarkerFaceColor', 'b');
% hold on
% loglog(xFit, yFit);
% ylim([1e-7, 1])
% xlabel('Average RMS WFE (\lambda_0)')
% ylabel('Average Null Depth')
% 
% % Format string for legend's "fit" item (print the fit equation as well)
% ftEq = sprintf('Y = %0.1f \\times X^{%0.1f}',10^P(end),P(1));
% % Add legend
% hleg = legend('Data',['Fit : ',ftEq], 'Location', 'northwest');


%-- Plot
figure('Color', 'White');
% Charge 2
errorbar(wfe_mn_ch2, nll_mn_ch2, nll_std_ch2, nll_std_ch2, wfe_std_ch2, wfe_std_ch2,'o', 'MarkerFaceColor', 'b');
hold on
if ~isApproxOnly
    plot(xFit_ch2, yFit_ch2, 'c');
end
if isFitANDApprox || isApproxOnly
    plot(xFit_ch2, yFitMan_ch2,'k'); % Manual fit
end
% Charge 1
errorbar(wfe_mn_ch1, nll_mn_ch1, nll_std_ch1, nll_std_ch1, wfe_std_ch1, wfe_std_ch1,'o', 'MarkerFaceColor', 'r');
if ~isApproxOnly
    plot(xFit_ch1, yFit_ch1,'m');
end
if isFitANDApprox || isApproxOnly
    plot(xFit_ch1, yFitMan_ch1,'g');  % Manual fit
end

%-- Format the plot
set(gca, 'XScale','log', 'YScale','log')
ylim([1e-7, 1])
xlabel('Average RMS WFE (\lambda_0)')
ylabel('Average Null Depth')
title('RMS WFE to Null')
grid on

% Format string for legend's "fit" item (print the fit equation as well)
ftEq_ch2 = sprintf('Y = %0.3f \\times X^{%0.2f}',10^P_ch2(end),P_ch2(1));
ftEq_ch1 = sprintf('Y = %0.3f \\times X^{%0.2f}',10^P_ch1(end),P_ch1(1));

manEq_ch2 = sprintf('Y = (%0.3f \\times X)^2',manSlp_ch2);
manEq_ch1 = sprintf('Y = (%0.3f \\times X)^2',manSlp_ch1);
% Add legend 
if isFitANDApprox
    % (fit AND manual approx.)
    hleg = legend('l=2 (1 STD err)',['Fit (l=2): ',ftEq_ch2],['Approx. (l=2): ',manEq_ch2],'l=1 (1 STD err)', ['Fit (l=1): ',ftEq_ch1],['Approx. (l=1): ',manEq_ch1], 'Location', 'northwest');
elseif isApproxOnly
    % (manual approx only)
    hleg = legend('l=2 (1 STD err)',['Approx. (l=2): ',manEq_ch2],'l=1 (1 STD err)', ['Approx. (l=1): ',manEq_ch1], 'Location', 'northwest');
else
    % (just the data and auto fits)
    hleg = legend('l=2 (1 STD err)',['Fit (l=2): ',ftEq_ch2],'l=1 (1 STD err)', ['Fit (l=1): ',ftEq_ch1], 'Location', 'northwest');
end

%% Zernike Analysis
%-- Load the data (using data from VFN server KPIC simulations)
% Average zernike values (by coeff)
zrn_mn = [ 0.004081, 
           1e-20;   %Actually 0 but set very small for log plotting
           1e-20;   %Actually 0 but set very small for log plotting
          -0.004112;
          -0.000163;
           0.002355;
           0.002326;
          -0.000780;
          -0.000125;
          -0.001437;
          -0.001171;
          -0.000745;
          -0.000539;
           0.001665;
           0.000910;
          -0.001893;
          -0.000059;
           0.000300;
          -0.000087;
           0.000069;
          -0.000487;
           0.002577;
          -0.000134;
          -0.000193;
          -0.000437;
          -0.000043;
          -0.000167;
          -0.000431;
          -0.000301; 
           0.001107;
           0.000168;
           0.000121;
          -0.000305;
           0.000498;
           0.000708;
           0.000629]';
       
% STD of zernike values (by coeff)
 zrn_std = [ 0.006391;
             0;
             0;
             0.009718;
             0.009878;
             0.010896;
             0.008137;
             0.008130;
             0.006505;
             0.008241;
             0.005947;
             0.007800;
             0.007251;
             0.006271;
             0.005998;
             0.004783;
             0.004910;
             0.005595;
             0.006083;
             0.005111;
             0.005128;
             0.005702;
             0.004465;
             0.004174;
             0.004651;
             0.004919;
             0.004202;
             0.005440;
             0.003936;
             0.003844;
             0.003712;
             0.003915;
             0.004461;
             0.004077;
             0.004694;
             0.004807]';

%-- Get indices of pertinent Zernikes
%- Convert from noll to radial/azimuthal indices:
% addpath('..\VFNlib'))
% for jj = 1:36
%     [n(jj),m(jj)] = zern_noll2index(jj);
% end
%- Extract m=1 (coma)  values
% iscoma = abs(m) == 1;
%- Extract m=2 (astig) values
% isastg = abs(m) == 2;
%- (Results from above without running that code)
iscoma = [0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0];
isastg = [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%-- Find all non-astig and non-coma terms
isrest = double(~(isastg | iscoma));
%- Convert 0s to nans for simple plotting
iscoma(iscoma==0) = nan;
isastg(isastg==0) = nan;
isrest(isrest==0) = nan;

%-- Plot
zrn_xax = 1:length(zrn_mn);
figure('Color', 'White');
% Errorbar version
% errorbar(zrn_xax, abs(zrn_mn.*iscoma), zrn_std.*iscoma, 'o', 'Color', 'r', 'MarkerFaceColor', 'r');
% hold on
% errorbar(zrn_xax, abs(zrn_mn.*isastg), zrn_std.*isastg, 'o', 'Color', 'b', 'MarkerFaceColor', 'b');
% errorbar(zrn_xax, abs(zrn_mn.*isrest), zrn_std.*isrest, 'o', 'Color', 'k');
% Scatter version
scatter(zrn_xax, abs(zrn_mn.*iscoma)/max(zrn_mn), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
hold on
scatter(zrn_xax, abs(zrn_mn.*isastg)/max(zrn_mn), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
scatter(zrn_xax, abs(zrn_mn.*isrest)/max(zrn_mn), 'o', 'MarkerEdgeColor', 'k');

%-- Format Plot
set(gca, 'YScale','log')
xlim([1 36])
ylim([1e-2 1])
xlabel('Noll Index')
ylabel('Average Zernike Coeff (peak norm)')
title('Zernike Coefficients')
legend('Coma Terms', 'Astig Terms', 'Other Terms')
grid on
grid minor
