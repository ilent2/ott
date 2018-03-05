function R = calc_rotation_matrix(  angle_vector )
% calc_rotation_matrix : finds coordinate transformation matrix
% for pure rotation.
%
% The rotation is specified by the angle vector, given by
% angle_vector = integral(0,t) angular_velocity_vector(t) dt
%
% x' = x * R; for row vector (or matrix of row vectors)
%
% NB: Since this  is an SI package, the angles are all in radians!
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

R = [ (cos(angle_vector(2))*cos(angle_vector(3))) ...
      (sin(angle_vector(3))) ...
      (-sin(angle_vector(2)));
   (-sin(angle_vector(3))) ...
      (cos(angle_vector(1))*cos(angle_vector(3))) ...
      (-sin(angle_vector(1)));
   (sin(angle_vector(2))) ...
      (sin(angle_vector(1))) ...
      cos(angle_vector(1))*cos(angle_vector(2)) ];   

return
