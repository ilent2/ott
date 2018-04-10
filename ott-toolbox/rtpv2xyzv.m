function [xv,yv,zv,x,y,z] = rtpv2xyzv(rv,thetav,phiv,r,theta,phi)
% rtpv2xyzv.m : convert a vector field of vector values at
%               points from spherical to cartesian
%
% Usage:
% [xv,yv,zv,x,y,z] = rtpv2xyzv(rv,thetav,phiv,r,theta,phi)
% or
% [vec_cart,pos_cart] = rtpv2xyzv(vec,pos)
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin < 6
   r = thetav(:,1);
   theta = thetav(:,2);
   phi = thetav(:,3);
   phiv = rv(:,3);
   thetav = rv(:,2);
   rv = rv(:,1);
end

% Convert points to cartesian coordinates
[x,y,z] = rtp2xyz(r,theta,phi);

%Calculate the Jacobian for spherical to cartesian unit vectors
%(transpose of catesian to spherical).
J=[sin(theta).*cos(phi),cos(theta).*cos(phi),-sin(phi);...
    sin(theta).*sin(phi),cos(theta).*sin(phi),cos(phi);...
    cos(theta),-sin(theta),zeros(size(theta))];

%Pack the spherical three vectors
rtpv=[rv,thetav,phiv];

%Separate the Jacobian and multiply for each unit vector.
xv = dot(J(1:length(theta),:),rtpv,2);
yv = dot(J(length(theta)+1:2*length(theta),:),rtpv,2);
zv = dot(J(2*length(theta)+1:3*length(theta),:),rtpv,2);

if nargout < 3
   xv = [ xv yv zv ];
   yv = [ x y z ];
end

return

