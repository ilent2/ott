function alpha = rotate_3x3tensor(ualpha, varargin)
% Apply a set of rotations to a 3x3 tensor.
%
% Applies the operation::
%
%   alpha = R*ualpha*inv(R)
%
% Where `R` are 3x3 rotation matrices and ualpha is a 3x3 tensor.
% Supports arrays of rotation matrices.
%
% Usage
%   alpha = rotate_3x3tensor(ualpha, R, ...)
%
%   alpha = rotate_3x3tensor(ualpha, 'direction', direction, ...)
%   computes the appropriate rotation matrix to rotate from the z-axis
%   to the 3xN matrix of directions.  This is useful for uniaxial materials.
%
% Parameters
%   - ualpha (3x3 numeric) -- Tensor to rotate.
%   - R (3x3N numeric) -- Rotation matrices.
%   - direction (3xN numeric) -- Cartesian vectors [x;y;z].
%   - sphdirection (2xN numeric) -- Spherical vectors [phi; theta]
%
% Optional named parameters
%   - inverse (logical) -- When true, returns the inverse tensor.
%     Default: false.

% Copyright 2018 Isaac Lenton (aka ilent2)
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.rotz;
import ott.utils.roty;

p = inputParser;
p.addOptional('rotation', []);
p.addParameter('direction', []);
p.addParameter('sphdirection', []);
p.addParameter('inverse', false);
p.parse(varargin{:});

assert(isnumeric(ualpha) && all(size(ualpha) == [3, 3]), ...
    'ualpha must be 3x3 matrix');

assert(~isempty(p.Results.rotation) + ~isempty(p.Results.direction) ...
   + ~isempty(p.Results.sphdirection) == 1, ...
    'Must supply one of direction, sphdirection or rotation');

if ~isempty(p.Results.direction)
  direction = p.Results.direction;
  assert(size(direction, 1) == 3, 'Direction must be 3xN matrix');

  % Calculate rotations off z-axis
  theta = atan2(direction(2, :), direction(1, :)) * 180/pi;
  phi = atan2(sqrt(direction(2, :).^2 + direction(1, :).^2), ...
      direction(3, :)) * 180/pi;

  % Calculate rotation matrices
  rotation = zeros(3, 3*size(phi, 2));
  for ii = 1:size(phi, 2)
    rotation(:, (1:3) + 3*(ii-1)) = rotz(theta(ii))*roty(phi(ii));
  end

elseif ~isempty(p.Results.sphdirection)
  assert(size(p.Results.sphdirection, 1) == 2, ...
      'sphdirection must be 2xN matrix');

  phi = p.Results.sphdirection(1, :);
  theta = p.Results.sphdirection(2, :);

  % Calculate rotation matrices
  rotation = zeros(3, 3*size(phi, 2));
  for ii = 1:size(phi, 2)
    rotation(:, (1:3) + 3*(ii-1)) = rotz(theta(ii))*roty(phi(ii));
  end

else
  rotation = p.Results.rotation;
  assert(size(rotation, 1) == 3 && mod(size(rotation, 2), 3) == 0, ...
      'Rotation must be 3x3N matrix');
end

% Calculate polarizabilities
alpha = zeros(size(rotation));
for ii = 1:size(alpha, 2)/3
  R = rotation(:, (1:3) + 3*(ii-1));
  alpha(:, (1:3) + 3*(ii-1)) = R*ualpha*R.';
end

% Calculate inverse if requested
if p.Results.inverse
  for ii = 1:size(alpha, 2)/3
    alpha(:, (1:3) + 3*(ii-1)) = inv(alpha(:, (1:3) + 3*(ii-1)));
  end
end

