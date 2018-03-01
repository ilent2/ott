function [M,N,M2,N2,M3,N3] = vswfsph2cart(n,m,kr,theta,phi,htype)
% vswfsph2cart.m : vector spherical harmonics
%                  spherical coordinate input, cartesian output
%
% Usage:
% [M,N] = vswfsph2cart(n,m,kr,theta,phi,type)
% or
% [M1,N1,M2,N2,M3,N3] = vswfsph2cart(n,m,kr,theta,phi)
%
% where
% kr, theta, phi are vectors of equal length, or scalar.
% type = 1 -> outgoing solution - h(1)
% type = 2 -> incoming solution - h(2)
% type = 3 -> regular solution - j (ie RgM, RgN)
%
% Scalar n,m for the moment.
% M,N are arrays of size length(vector_input) x 3
%
% The three components of each input vector are [kr,theta,phi]
% The three components of each output vector are [x,y,z]
%
% "Out of range" n and m result in return of [0 0 0]
%
% PACKAGE INFO


if nargin < 6
   htype = 0;
end

[M,N,M2,N2,M3,N3] = vswf(n,m,kr,theta,phi,htype);

% Convert to cartesian coordinates
[x,y,z] = rtp2xyz(kr,theta,phi);
theta_hat_x = cos(theta) .* cos(phi);
theta_hat_y = cos(theta) .* sin(phi);
theta_hat_z = -sin(theta);
phi_hat_x = -sin(phi);
phi_hat_y = cos(phi);
phi_hat_z = 0;
r_hat_x = x./kr;
r_hat_y = y./kr;
r_hat_z = z./kr;

if length(M) > 1
Mr = M(:,1);
Mtheta = M(:,2);
Mphi = M(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
M = [ Mx My Mz ];
end

if length(N) > 1
Mr = N(:,1);
Mtheta = N(:,2);
Mphi = N(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
N = [ Mx My Mz ];
end

if length(M2) > 1
Mr = M2(:,1);
Mtheta = M2(:,2);
Mphi = M2(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
M2 = [ Mx My Mz ];
end

if length(N2) > 1
Mr = N2(:,1);
Mtheta = N2(:,2);
Mphi = N2(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
N2 = [ Mx My Mz ];
end

if length(M3) > 1
Mr = M3(:,1);
Mtheta = M3(:,2);
Mphi = M3(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
M3 = [ Mx My Mz ];
end

if length(N3) > 1
Mr = N3(:,1);
Mtheta = N3(:,2);
Mphi = N3(:,3);
Mx = Mr .* r_hat_x + Mtheta .* theta_hat_x + Mphi .* phi_hat_x;
My = Mr .* r_hat_y + Mtheta .* theta_hat_y + Mphi .* phi_hat_y;
Mz = Mr .* r_hat_z + Mtheta .* theta_hat_z + Mphi .* phi_hat_z;
N3 = [ Mx My Mz ];
end

return
