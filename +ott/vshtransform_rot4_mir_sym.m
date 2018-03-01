function [b,c] = vshtransform_rot4_mir_sym(f,n,m,nmax,theta,phi);
% vshtransform_rot4_mir_sym.m : vector spherical harmonic transform
% Redundant modes due to symmetry are not calculated 
%
% Finds coefficients for the expansion of a 2D angular
% function with a 3D vector value in terms spherical harmonics.
%
% The expansion coefficients a(n,m) are defined by
% f = sum( a(n,m) * P(n,m) + b(n,m) * B(n,m) + c(n,m) * C(n,m) )
%
% Usage:
% a = vshtransform(f,nmax,theta,phi)
% where each row of f is the three vector components
% in spherical coordinates
%
% The theta and phi vectors really should be produced using
% [theta,phi] = angulargrid(ntheta,nphi);
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

% The coefficients are found using the formula
%
% a(n) = integral f(theta,phi).P*(n,m,theta,phi) dOmega
%        / integral P(n,m,theta,phi).P*(n.m.theta,phi) dOmega
%
% Note that the denominator = 1 (for P) and n(n+1) (for B,C).

total_harmonics = nmax * ( nmax + 2 ) ;

%a = zeros(total_harmonics,1); not required
b = zeros(total_harmonics,1);
c = zeros(total_harmonics,1);

ntheta = max(find(phi==phi(1)));
nphi = length(theta)/ntheta;
[thetav,phiv] = angulargrid(ntheta,nphi,1); 
thetam = reshape(theta,ntheta,nphi);
phim = reshape(phi,ntheta,nphi);
sint = sin(thetav(:));

[nn, mm, cii] = nm_oct(n,m,nmax);
kk =0;
for ii = 1:length(nn)
  ci = cii(ii);
    
  [B,C,P] = vsh(nn(ii),mm(ii),theta,phi);

  dotf = dot(B,f,2);
  dotf = reshape(dotf,ntheta,nphi);
  phi_int = trapz(phiv,dotf,2);
  surface_integral = trapz(thetav,phi_int.*sint);
  b(ci) = surface_integral/(n*(n+1));

  dotf = dot(C,f,2);
  dotf = reshape(dotf,ntheta,nphi);
  phi_int = trapz(phiv,dotf,2);
  surface_integral = trapz(thetav,phi_int.*sint);
  c(ci) = surface_integral/(n*(n+1));
end

return
