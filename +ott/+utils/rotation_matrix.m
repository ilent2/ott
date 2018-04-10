function R = rotation_matrix(rot_axis,rot_angle)
% ROTATION_MATRIX calculates rotation matrix using Euler-Rodrigues formula.
%
% R = rotation_matrix( axis, angle ) calculates the rotation about
% axis by angle (in radians).
%
% R = rotation_matrix( axis_angle ) calculates the rotation about vector
% axis_angle, the angle is specified as the length of the vector (in radians).
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin < 2
	rot_angle= norm(rot_axis);
    	rot_axis = rot_axis / norm(rot_axis);
end

if rot_angle==0 || norm(rot_axis)==0
	R=eye(3);
	return
end

rot_axis=rot_axis / norm(rot_axis);

% Euler-Rodrigues formula, might be unstable for small angles :(

ux = [ 0 -rot_axis(3) rot_axis(2);
    rot_axis(3) 0 -rot_axis(1);
    -rot_axis(2) rot_axis(1) 0 ];

R = eye(3) + sin(rot_angle) * ux + ( 1 - cos(rot_angle) ) * ux * ux;
