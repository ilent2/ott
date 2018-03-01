function [x,y,z] = rtp2xyz(r,theta,phi)
% rtp2xyz.m
% Coordinate transformation from spherical to cartesian
% theta is the polar angle, measured from the +z axis,
% and varies from 0 to pi
% phi is the azimuthal angle, measured from the +x axis, increasing
% towards the +y axis, varying from 0 to 2*pi
%
% Usage:
% [x,y,z] = rtp2xyz(r,theta,phi);
% where x, y, z, r, theta, phi are all scalars or
% equal-length vectors
% or
% x = rtp2xyz(r);
% where x = [ x y z ] and r = [ r theta phi ]
%
% Angles are in radians
%
% PACKAGE INFO

if nargin == 1
   theta = r(:,2);
   phi = r(:,3);
   r = r(:,1);
end

z = r .* cos(theta);
xy = r .* sin(theta);

x = xy .* cos(phi);
y = xy .* sin(phi);

if nargout == 1
   x = x(:);
   y = y(:);
   z = z(:);
   x = [ x y z ];
end

return
