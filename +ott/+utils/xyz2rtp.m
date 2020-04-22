function varargout = xyz2rtp(x,y,z)
% Coordinate transformation from Cartesian to Spherical
%
% Usage
%   [r, theta, phi] = xyz2rtp(x, y, z)
%
%   rtp = xyz2rtp(xyz) as above but using 3xN matrices for input/output.
%
% Parameters
%   - x,y,z -- Cartesian coordinates
%   - r     -- radial distance [0, Inf)
%   - theta -- polar angle, measured from +z axis [0, pi]
%   - phi   -- azimuthal angle, measured from +x towards +y axes [0, 2*pi)

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 1
  assert(isnumeric(x) && size(x, 1) == 3, 'xyz must be 3xN numeric matrix');

  y = x(2, :).';
  z = x(3, :).';
  x = x(1, :).';

elseif nargin == 3
  % Nothing to do
else
  error('Must supply either xyz or x, y, z');
end

xy = sqrt( x.*x + y.*y );
theta = mod(atan2(xy,z)+2*pi,2*pi);
phi = mod(atan2(y,x)+2*pi,2*pi);
r = sqrt( x.*x + y.*y + z.*z );

if nargout == 1
  varargout{1} = [r(:), theta(:), phi(:)].';
else
  varargout{1} = r;
  varargout{2} = theta;
  varargout{3} = phi;
end

