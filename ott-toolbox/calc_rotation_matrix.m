function R = calc_rotation_matrix(  angle_vector )
%CALC_ROTATION_MATRIX finds rotation matrix
%
% R = CALC_ROTATION_MATRIX(angle_vector) calculates rotation matrix.
% The rotation is specified by the angle vector, given by
% angle_vector = integral(0,t) angular_velocity_vector(t) dt
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:calc_rotation_matrix:depreciated', ...
    ['calc_rotation_matrix.m will be depreciated in ott1.4. ' ...
     'rotation_matrix.m can be used instead.']);

ott_warning('internal');

if length(angle_vector)==2
    angle_vector=[0;angle_vector(:)];
end

R = [ (cos(angle_vector(2))*cos(angle_vector(3))) ...
      (sin(angle_vector(3))) ...
      (-sin(angle_vector(2)));
   (-sin(angle_vector(3))) ...
      (cos(angle_vector(1))*cos(angle_vector(3))) ...
      (-sin(angle_vector(1)));
   (sin(angle_vector(2))) ...
      (sin(angle_vector(1))) ...
      cos(angle_vector(1))*cos(angle_vector(2)) ];   

ott_warning('external');

return
