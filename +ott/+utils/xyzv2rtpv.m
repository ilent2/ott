function [rv,thetav,phiv,r,theta,phi] = xyzv2rtpv(xv,yv,zv,x,y,z)
% Cartesian to spherical vector field conversion
%
% Usage
%   [rv,thetav,phiv,r,theta,phi] = XYZV2RTPV(xv,yv,zv,x,y,z)
%
%   [vec_sph,pos_sph] = XYZV2RTPV(vec_cart, pos_cart)
%   As above, but with 3xN matrices as inputs and outputs.
%
% Parameters
%   - r      -- radial distance [0, Inf)
%   - theta  -- polar angle, measured from +z axis [0, pi]
%   - phi    -- azimuthal angle, measured from +x towards +y axes [0, 2*pi)
%   - rv     -- Vector radial distance
%   - thetav -- Vector polar angle
%   - phiv   -- Vector azimuthal angle
%   - x,y,z  -- Cartesian position coordinates
%   - xv,yv,zv -- Cartesian vector coordinates
%
% See also rtpv2xyzv and xyz2rtp.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 2

  assert(size(xv, 1) == 3, 'pos must be 3xN matrix');
  assert(size(yv, 1) == 3, 'vec must be 3xN matrix');
  assert(size(xv, 2) == size(yv, 2), ...
      'Number of points in vec and pos must match');

   x = yv(1,:).';
   y = yv(2,:).';
   z = yv(3,:).';
   zv = xv(3,:).';
   yv = xv(2,:).';
   xv = xv(1,:).';

elseif nargin == 6

  % Ensure inputs are columns
  xv = xv(:);
  yv = yv(:);
  zv = zv(:);
  x = x(:);
  y = y(:);
  z = z(:);

  % Check size
  N = unique([numel(xv), numel(yv), numel(zv), ...
      numel(x), numel(y), numel(z)]);
  assert(numel(N) == 1, 'Inputs must all be same length');

else
  error('Must supply either 2 or 6 inputs');
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
   rv = [ rv(:) thetav(:) phiv(:) ].';
   thetav = [ r(:) theta(:) phi(:) ].';
end

