function [nn,mm,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
% Finds VSWF representation of plane wave
% [n,m,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
%
% theta and phi (in radians) give the direction of propagation
% of the plane wave. +z direction is theta = 0, phi = any value
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

E = [ 0 Etheta Ephi ];

mm = [];
nn = [];
a = [];
b = [];

for n = 1:nmax
  
   Nn = sqrt(1/(n*(n+1)));
  
   for m = -n:n
      
      [B,C,P] = vsh(n,m,theta,phi);
      
      % dot product function takes complex conjugate of first term
      newa = 4*pi * i^n * Nn * dot(C,E);
      newb = 4*pi * i^(n-1) * Nn * dot(B,E);
      
      % Get rid of any values too close to zero
      if abs(newa) + abs(newb) > 1e-14
         nn = [ nn n ];
         mm = [ mm m ];
         a = [ a newa ];
         b = [ b newb ];
      end
      %fprintf(1,'Done n = %d, m = %d, a = %g+%gi, b = %g+%gi.\n',...
      %   n,m,real(newa),imag(newa),real(newb),imag(newb));
      
   end
   
end

return
