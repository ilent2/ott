function R = rotation_matrix(rot_axis,rot_angle)
% ROTATION_MATRIX calculates rotation matrix using Euler-Rodrigues formula.
%
% R = rotation_matrix( axis, angle ) calculates the rotation about
% axis by angle (in radians).
%
% R = rotation_matrix( axis_angle ) calculates the rotation about vector
% axis_angle, the angle is specified as the length of the vector (in radians).

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin < 2
	rot_angle = vecnorm(rot_axis);
  rot_axis = rot_axis ./ rot_angle;
end

assert(isnumeric(rot_angle) && isscalar(rot_angle), ...
    'rot_angle must be numeric scalar');
assert(isnumeric(rot_axis) && isvector(rot_axis) ...
  && numel(rot_axis) == 3, ...
  'rot_axis should be 3 element numeric vector');

% Check for small angle
tol = 1.0e-8;
if abs(rot_angle) < tol
  R = eye(3);
  return;
end

% Check axis is unit vector
if abs(vecnorm(rot_axis)-1) > 2*eps(1)
  warning('ott:utils:rotation_matrix:unit_vector', ...
    'Rotation axis (with two inputs) should be unit vector, rescaling');
  rot_axis = rot_axis ./ vecnorm(rot_axis);
end

% Euler-Rodrigues formula, might be unstable for small angles :(

ux = [ 0 -rot_axis(3) rot_axis(2);
    rot_axis(3) 0 -rot_axis(1);
    -rot_axis(2) rot_axis(1) 0 ];

R = eye(3) + sin(rot_angle) * ux + ( 1 - cos(rot_angle) ) * ux * ux;
