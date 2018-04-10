function R = rotation_matrix(rot_axis,rot_angle)
% rotation_matrix.m
%
% Usage:
% R = rotation_matrix( axis, angle );
% R = rotation_matrix( axis_angle );
% where
% axis is a vector giving the direction of the rotation axis, any length,
% angle is the angle of rotation, in radians,
% axis_angle is a vector giving the rotation axis, of magnitude equal to
% the angle of rotation.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

if nargin < 2
	rot_angle= norm(rot_axis);
    	rot_axis = rot_axis / norm(rot_axis);
end

if rot_angle==0|norm(rot_axis)==0
	R=eye(3);
	return
end

rot_axis=rot_axis / norm(rot_axis);

% Euler-Rodrigues formula, might be unstable for small angles :(

ux = [ 0 -rot_axis(3) rot_axis(2);
    rot_axis(3) 0 -rot_axis(1);
    -rot_axis(2) rot_axis(1) 0 ];

R = eye(3) + sin(rot_angle) * ux + ( 1 - cos(rot_angle) ) * ux * ux;

return
