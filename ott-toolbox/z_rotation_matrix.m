function [Ryz]=z_rotation_matrix(theta,phi,tolerance)
% z_rotation_matrix.m: Generates the rotation matrix for use in the wigner
%                       rotation function using the product of Ry, Rz.
% 
% USAGE:
%
% [Ryz]=z_rotation_matrix(theta,phi,tolerance)
% 
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
% PACKAGE INFO

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
