function [r,theta,phi] = xyz2rtp(x,y,z)
% xyz2rtp.m
% Coordinate transformation from cartesian to spherical
% theta is the polar angle, measured from the +z axis,
% and varies from 0 to pi
% phi is the azimuthal angle, measured from the +x axis, increasing
% towards the +y axis, varying from 0 to 2*pi
%
% Usage:
% [r,theta,phi] = xyz2rtp(x,y,z);
% where x, y, z, r, theta, phi are all scalars or
% equal length vectors
% or
% r = xyz2rtp(x);
% where x = [ x y z ] and r = [ r theta phi ]
%
% Angles are in radians
%
% PACKAGE INFO

if nargin == 1
   y = x(:,2);
   z = x(:,3);
   x = x(:,1);
end

xy = sqrt( x.*x + y.*y );
theta = mod(atan2(xy,z)+2*pi,2*pi);
phi = mod(atan2(y,x)+2*pi,2*pi);
r = sqrt( x.*x + y.*y + z.*z );

if nargout == 1
   r = r(:);
   theta = theta(:);
   phi = phi(:);
   r = [ r theta phi ];
end

return
