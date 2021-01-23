function [ OUT ] = getWedgeOffset(theta, lambda)
%Returns an angular offset for the phase of the pupil plane
%induced by the addition of an optical wedge with index of refraction (n)
%and angle between surfaces (theta)

n = getRefractiveIndex('NBK7', lambda);

[OUT] = (n-1)*theta;

