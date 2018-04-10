function [Ryz]=z_rotation_matrix(theta,phi,tolerance)
% Z_ROTATION_MATRIX generates the rotation matrix for use in the wigner
% rotation function using the product of Ry, Rz.
% 
% [Ryz]=Z_ROTATION_MATRIX(theta,phi,tolerance)
% phi is the azimuthal angle from the +x.
% theta is the elevation angle from +z.
% 
% Can include a tolerance so that round off zeros are zero using the most
% significant value in the matrix as the scaling value.
%
% WARNING:  Trancation with small perturbation angles is not recommended.
% This is not a general rotation matrix, it only works if you start on the
% z-axis.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

%warning('ott:z_rotation_matrix:move','this function will move to ott.utils.z_rotation_matrix');

warning('ott:calc_rotation_matrix:depreciated', ...
    ['calc_rotation_matrix.m will be depreciated in ott1.4. ' ...
     'rotation_matrix.m can be used instead.']);
     
ct=cos(theta);
st=sin(theta);
cp=cos(phi);
sp=sin(phi);
% % Ryz = Rz*Ry
Ryz =   [ct*cp,   -sp,    st*cp; ...
        ct*sp,    cp,     st*sp; ...
        -st,      0,      ct];
    
if nargin>2
    Ryz(abs(Ryz)<max(abs(Ryz(:)))*tolerance)=0;
end
