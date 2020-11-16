function rtp = paraxial2rtp(xy, mapping, direction)
% Convert from paraxial coordinates to spherical coordinates.
%
% Coordinates have the range ``r > 0; 0 <= t <= pi; 0 <= p < 2*pi``.
% Values outside the paraxial hemisphere are set to NaN.
%
% Usage
%   rtp = paraxial2rtp(xy, mapping, direction)
%
% Parameters
%   - xy (2xN numeric) -- Paraxial coordinates
%
%   - mapping (enum) -- Mapping from theta-phi to far-field.
%     Must be one of 'sin', 'tan' or 'theta'.
%
%   - direction (enum | 2 numeric) -- Mapping direction.
%     Can be enum: ``'pos'`` or ``'neg'`` for +z and -z hemisphere, or
%     `[theta, phi]` specifying normal direction.
%
% See also :func:`rtp2paraxial` and :func:`xyz2rtp`.

% Copyright 2020 Isaac Lenton (aka ilent2)
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

phi = atan2(xy(2, :), xy(1, :))+pi;
rr = vecnorm(xy, 2, 1);
switch mapping
  case 'sin'
    theta = asin(rr);
  case 'tan'
    theta = atan(rr);
  case 'theta'
    theta = rr;
  otherwise
    error('Unknown mapping argument value, must be sin, tan or theta');
end

% Remove imaginary theta/phi
phi(imag(theta) ~= 0) = nan;
theta(imag(theta) ~= 0) = nan;

% Adjust coordinates for direction
if ischar(direction)
  if strcmpi(direction, 'neg')
    theta = pi - theta;
  elseif strcmpi(direction, 'pos')
    % Nothing to do
  else
    error('ott:utils:paraxial2rtp:unknown_direction_enum', ...
        'Unknown direction enum specified');
  end
elseif numel(direction) == 2
  theta = theta + direction(1);
  phi = phi + direction(2);
end

rtp = [zeros(size(theta)); theta; phi];

