function [x,y,z] = rtp2xyz(r,theta,phi)
% RTP2XYZ coordinate transformation from spherical to cartesian
% theta is the polar angle, measured from the +z axis,
% and varies from 0 to pi
% phi is the azimuthal angle, measured from the +x axis, increasing
% towards the +y axis, varying from 0 to 2*pi
%
% [x,y,z] = RTP2XYZ(r,theta,phi) takes vectors or scalars for:
%   r      radial distance [0, Inf)
%   theta  polar angle, measured from +z axis [0, pi]
%   phi    azimuthal angle, measured from +x axis [0, 2*pi)
%
% [x,y,z] = RTP2XYZ(r) same as above but with the coordinate
% packed into the vector/matrix r = [ r theta phi ].
%
% x = RTP2XYZ(...) same as above with the result packed into
% the vector/matrix x = [ x y z ].
%
% PACKAGE INFO

warning('this function will move to ott.utils.rtp2xyz');

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
