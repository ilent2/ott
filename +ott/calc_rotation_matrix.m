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
% PACKAGE INFO

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

return
