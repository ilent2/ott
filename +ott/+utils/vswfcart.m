function varargout = vswfcart(n,m,kr,theta,phi,htype)
% VSWFCART vector spherical harmonics spherical coordinate input,
% cartesian output.
%
% [M1,N1,M2,N2,M3,N3] = VSWFCART(n,m,kr,theta,phi) calculates the
% outgoing M1,N1, incomming M2,N2 and regular M3,N3 VSWF.
% kr, theta, phi are vectors of equal length, or scalar.
%
% [M,N] = VSWFCART(n,m,kr,theta,phi,type) calculates only the
% requested VSWF, where type is
%     1 -> outgoing solution - h(1)
%     2 -> incoming solution - h(2)
%     3 -> regular solution - j (ie RgM, RgN)
%
% Scalar n,m for the moment.
% M,N are arrays of size length(vector_input) x 3
%
% The three components of each input vector are [kr,theta,phi]
% The three components of each output vector are [x,y,z]
%
% "Out of range" n and m result in return of [0 0 0]
%
% At the coordinate origin (kr == 0) we use only theta/phi.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*
ott.warning('internal');

if nargin < 6
   htype = 0;
end

[varargout{1:nargout}] = vswf(n,m,kr,theta,phi,htype);

% Convert to cartesian coordinates
[x,y,z] = rtp2xyz(kr,theta,phi);
theta_hat_x = cos(theta) .* cos(phi);
theta_hat_y = cos(theta) .* sin(phi);
theta_hat_z = -sin(theta);
phi_hat_x = -sin(phi);
phi_hat_y = cos(phi);
phi_hat_z = 0;

% This avoids NaNs and makes debugging a little easier
% May also provide performance improvements (untested speculation)
kr_safe = kr;
kr_safe(kr == 0) = 1.0;

% TODO: Discuss if we actually need or even should use this
r_hat_x = x./kr_safe;
r_hat_y = y./kr_safe;
r_hat_z = z./kr_safe;

% The following should be equivilant to
%   M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
%     cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
%     -sin(phi) cos(phi) 0].';
%
% This requires modifying at least r_hat_z when kr == 0.0
r_hat_x(kr == 0) = sin(theta(kr == 0)).*cos(phi(kr == 0));
r_hat_y(kr == 0) = sin(theta(kr == 0)).*sin(phi(kr == 0));
r_hat_z(kr == 0) = cos(theta(kr == 0));

for ii = 1:nargout
  
  if isempty(varargout{ii})
    continue;
  end
  
  % utils.vswf might return multiple points, split acordingly
  npts = size(varargout{ii}, 2)/3;
  Mr = varargout{ii}(:, 1:npts);
  Mtheta = varargout{ii}(:, npts+(1:npts));
  Mphi = varargout{ii}(:, 2*npts+(1:npts));
  
  Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
  My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
  Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
  
  varargout{ii} = [Mx My Mz];
  
end

ott.warning('external');
