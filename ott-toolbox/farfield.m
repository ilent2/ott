function [E,H] = farfield(n,m,a,b,p,q,theta,phi)
% farfield.m
% Finds far field given VSWF expansion coefficients
%
% usage
% E = farfield(n,m,a,b,p,q,theta,phi)
%
% each row of E is the field (in spherical coordinates) in the
% (theta,phi) direction (assuming a distance scaling factor of kr)
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

[theta,phi] = matchsize(theta,phi);
E = zeros(length(theta),3);
H = zeros(length(theta),3);

a = threewide(a);
b = threewide(b);
p = threewide(p);
q = threewide(q);

for nn = 1:length(n)
   
   [B,C,P] = vsh(n(nn),m(nn),theta,phi);
   Nn = 1/sqrt(n(nn)*(n(nn)+1));
   
   E = E + Nn * (i)^(n(nn)+1) * a(nn).*C + Nn * (i)^n(nn) * b(nn).*B + ...
      Nn * (-i)^(n(nn)+1) * p(nn).*C + Nn * (-i)^n(nn) * q(nn).*B;
%   H = H + a.*N2 + b.*M2 + p.*N1 + q.*M1;
   
end

% SI-ify units of H
H = H * -1i;

return
