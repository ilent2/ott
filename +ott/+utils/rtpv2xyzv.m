function [xv,yv,zv,x,y,z] = rtpv2xyzv(rv,thetav,phiv,r,theta,phi)
% Spherical to Cartesian vector field conversion
%
% Usage
%   [xv,yv,zv,x,y,z] = RTPV2XYZV(rv,thetav,phiv,r,theta,phi)
%
%   [vec_cart,pos_cart] = rtpv2xyzv(vec,pos)
%   As above, but with 3xN matrices as inputs and outputs.
%   vec can also be 3xNxM, but pos must be a 3xN matrix.
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
% See also rtp2xyz and xyzv2rtpv.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 2
  
  assert(size(thetav, 1) == 3 && ismatrix(thetav), 'pos must be 3xN matrix');
  assert(size(rv, 1) == 3, 'vec must be 3xNxM matrix');
  assert(size(rv, 2) == size(thetav, 2), ...
      'Second dimension of vec and pos must match');
  
   r = thetav(1,:).';
   theta = thetav(2,:).';
   phi = thetav(3,:).';
   
   szrv = size(rv);
   phiv = reshape(rv(3,:), [numel(r), 1, szrv(3:end)]);
   thetav = reshape(rv(2,:), [numel(r), 1, szrv(3:end)]);
   rv = reshape(rv(1,:), [numel(r), 1, szrv(3:end)]);
elseif nargin == 6

  % Ensure inputs are columns
  rv = rv(:);
  thetav = thetav(:);
  phiv = phiv(:);
  r = r(:);
  theta = theta(:);
  phi = phi(:);

  % Check size
  N = unique([numel(rv), numel(thetav), numel(phiv), ...
      numel(r), numel(theta), numel(phi)]);
  assert(numel(N) == 1, 'Inputs must all be same length');
else
  error('Must supply 2 or 6 arguments as input');
end

% Convert points to cartesian coordinates
[x,y,z] = ott.utils.rtp2xyz(r,theta,phi);

%Calculate the Jacobian for spherical to cartesian unit vectors
%(transpose of catesian to spherical).
J=[sin(theta).*cos(phi),cos(theta).*cos(phi),-sin(phi);...
    sin(theta).*sin(phi),cos(theta).*sin(phi),cos(phi);...
    cos(theta),-sin(theta),zeros(size(theta))];

% Pack the spherical three vectors
rtpv = cat(2, rv, thetav, phiv);

% Separate the Jacobian and multiply for each unit vector.
xv = sum(J(1:length(theta),:) .* rtpv, 2);
yv = sum(J(length(theta)+1:2*length(theta),:) .* rtpv, 2);
zv = sum(J(2*length(theta)+1:3*length(theta),:) .* rtpv, 2);

if nargout < 3
   xv = permute(cat(2, xv, yv, zv), [2, 1 3:ndims(rtpv)]);
   yv = [ x(:), y(:), z(:) ].';
end

