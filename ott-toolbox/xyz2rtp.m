function [r,theta,phi] = xyz2rtp(x,y,z)
% XYZ2RTP coordinate transformation from cartesian to spherical
%   r      radial distance [0, Inf)
%   theta  polar angle, measured from +z axis [0, pi]
%   phi    azimuthal angle, measured from +x towards +y axes [0, 2*pi)
%
% [r,theta,phi] = XYZ2RTP(x,y,z) takes vectors or scalars outputs
% the spherical coordinates as vectors/scalars of the same size.
%
% [r,theta,phi] = XYZ2RTP(x) same as above but with the coordinate
% packed into the vector/matrix x = [ x y z ].
%
% r = XYZ2RTP(...) same as above with the result packed into
% the vector/matrix r = [ r theta phi ].
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

%warning('ott:xyz2rtp:move','this function will move to ott.utils.xyz2rtp');

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
