function [radvg, rvec] = radialAverage(img, cent, radpts, angpts)
% VFN_An_radAverage Return vector with average radial profile of an image
%
%   - This uses polarTransform() to tranform to polar coordinates then 
%       takes the average azimuthally
%   - This function automatically takes the average out to the largest 
%       radius possible without clipping the edge of the image.
%
%   ** NOTE: due to the way polarTransform centers the transformation, your
%       image should have odd length in all dimensions with the center value at
%       the middle. (ex: 21x21 w/ cent=[11,11] instead of 20x20 w/ cent=[10,10])
%       This is not a problem for VFN data since it is output in this way
%       automatically by the scanning script.
%   
%   [radvg, rvec] = VFN_An_radAverage(img, cent, radpts, angpts)
%     Azimuthally average s.t. mean radial profile is provided.
%     - 'img' image to average. Should be odd sized (see note above).
%     - 'cent' OPTIONAL: indices for origin of polar coordinates
%               Should be: [rowindex, colindex]
%               If not provided, cent = the center of the image
%     - 'radpts' OPTIONAL: number of radial points for interpolation.
%               If not provided, radpts = number of points between center 
%               and edge of image.
%     - 'angpts' OPTIONAL: number of azimuthal points for interpolation.
%               If not provided, angpts = 360
%
%     Returns:
%     - 'radvg' average radial profile (vector)
%     - 'rvec'  vector of radial positions in pixels
%
%   Examples:
%      [radvg, rvec] = VFN_An_radAverage(img)
%         Returns a vector containing the average radial profile of img 
%         centered at the center of the image. Also returns the vector of 
%         radial positions.
% 
%      [radvg, rvec] = VFN_An_radAverage(img, [22, 20])
%         Returns a vector containing the average radial profile of img 
%         centered at (row=22,col=20). Also returns the vector of radial 
%         positions.
%
%      [radvg, rvec] = VFN_An_radAverage(img, [22, 20], 100)
%         Same as above but with 100 interpolated points in radial profile
%
%      [radvg, rvec] = VFN_An_radAverage(img, [22, 20], 100)
%         Same as above but with 360 interpolated points in azimuthal dir.

if nargin <= 1
    % Center not provided: Get center of matrix (ceil odd-sized matrices)
    cent = ceil(size(img)/2);
end

% Get maximum radius without clipping on edge of image (so that average is 
%   always of the same number of points)
    % max(cent) takes the center value farthest from an edge which is then 
    % subtracted from both edges to see which is closest. Then take the min
    % between this value and the the cent to see what is closest to an edge
rmax = min([size(img)-max(cent), cent]);

if nargin <= 2
    % radpts not provided: use default which is rmax
    radpts = rmax;
end

if nargin <= 3
    % angpts not provided: use default which is 360
    angpts = 360;
end

% Perform polar transform; omit vector 
[pts, rvec, ~] = polarTransform(img, cent, rmax, radpts, angpts, 'linear');

% Calculate azimuthal average to return average radial profile
radvg = mean(pts,2);
end