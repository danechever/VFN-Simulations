function [ OUT ] = getWedgeTilt(wedge_mat, wedge_angle, lambda)
%Returns an angular offset for the phase of the pupil plane
%induced by the addition of an optical wedge with index of refraction (n)
%and angle between surfaces (theta)

n = getRefractiveIndex(wedge_mat, lambda);

[OUT] = (n-1)*wedge_angle;

