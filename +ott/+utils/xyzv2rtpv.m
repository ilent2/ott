function [rv,thetav,phiv,r,theta,phi] = xyzv2rtpv(xv,yv,zv,x,y,z)
% XYZV2RTPV cartiesian to spherical vector field conversion
%
% [rv,thetav,phiv,r,theta,phi] = XYZV2RTPV(xv,yv,zv,x,y,z)
%
% [vec_sph,pos_sph] = XYZV2RTPV(vec_cart,pos_cart)
%
% See also rtpv2xyzv and xyz2rtp.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott.warning('internal');

if nargin < 6
   x = yv(:,1);
   y = yv(:,2);
   z = yv(:,3);
   zv = xv(:,3);
   yv = xv(:,2);
   xv = xv(:,1);
end

% Convert points to spherical coordinates
[r,theta,phi] = ott.utils.xyz2rtp(x,y,z);

%Jacobian for Cartesian to spherical unit vectors
J=[sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta);...
    cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta);...
    -sin(phi),cos(phi),zeros(size(theta))];

%Pack Cartesian vector field
xyzv=[xv,yv,zv];

%Separate the Jacobian and multiply for each unit vector.
rv = dot(J(1:length(theta),:),xyzv,2);
thetav = dot(J(length(theta)+1:2*length(theta),:),xyzv,2);
phiv = dot(J(2*length(theta)+1:3*length(theta),:),xyzv,2);

if nargout < 3
   rv = [ rv thetav phiv ];
   thetav = [ r theta phi ];
end

ott.warning('external');
