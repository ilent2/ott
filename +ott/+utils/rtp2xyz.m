function varargout = rtp2xyz(r, theta, phi)
% Coordinate transformation from spherical to Cartesian
%
% Usage
%   [x, y, z] = rtp2xyz(r, theta, phi)
%
%   xyz = rtp2xyz(rtp) as above but using 3xN matrices for input/output.
%
% Parameters
%   - r     -- radial distance [0, Inf)
%   - theta -- polar angle, measured from +z axis [0, pi]
%   - phi   -- azimuthal angle, measured from +x towards +y axes [0, 2*pi)
%   - x,y,z -- Cartesian coordinates

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin == 1
  assert(isnumeric(r) && size(r, 1) == 3, 'rtp must be 3xN numeric matrix');

  theta = r(2, :).';
  phi = r(3, :).';
  r = r(1, :).';

elseif nargin == 3
  % Nothing to do
else
  error('Must supply either rtp or r, theta, phi');
end

z = r .* cos(theta);
xy = r .* sin(theta);

x = xy .* cos(phi);
y = xy .* sin(phi);

if nargout == 1
  varargout{1} = [x(:), y(:), z(:)].';
else
  varargout{1} = x;
  varargout{2} = y;
  varargout{3} = z;
end

