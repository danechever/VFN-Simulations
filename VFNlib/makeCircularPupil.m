function [ PUPIL ] = makeCircularPupil( rad, Ngrid )
%makeCircularPupil Summary of this function goes here
%   Detailed explanation goes here

    [X,Y] = meshgrid(-Ngrid/2:Ngrid/2-1);
    [~,RHO] = cart2pol(X,Y);

    PUPIL = exp(-(RHO/rad).^2000);
end

