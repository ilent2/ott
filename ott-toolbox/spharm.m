function [Y,Ytheta,Yphi] = spharm(n,m,theta,phi)
% spharm.m : scalar spherical harmonics and
%            angular partial derivatives
%
% Usage:
% Y = spharm(n,m,theta,phi)
% or
% [Y,dY/dtheta,1/sin(theta)*dY/dphi] = 
%
% Scalar n,m for the moment.
% Y is a vector of length(theta,phi)
%
% "Out of range" n and m result in return of Y = 0
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

[theta,phi] = matchsize(theta,phi);
input_length = length(theta);

if abs(m) > n | n < 0
   Y = zeros(input_length,1);
   Ytheta = zeros(input_length,1);
   Yphi = zeros(input_length,1);
   return
end

pnm = legendre(n,cos(theta),'sch');
pnm = pnm(abs(m)+1,:).';

% Why is this needed? Better do it, or m = 0 square integrals
% are equal to 1/2, not 1.
% This is either a bug in legendre or a mistake in the docs for it!
% Check this if MATLAB version changes! (Version 5.X)
if m == 0
   pnm = pnm * sqrt(2);
end

if m < 0
   pnm = (-1)^m * pnm;
end

expphi = exp(i*m*phi);

N = sqrt((2*n+1)/(8*pi));

Y = N * pnm .* expphi;

% Do we want to calculate the derivatives?
if nargout <= 1
   % Doesn't look like it
   return
end

% We use recursion relations to find the derivatives, choosing
% ones that don't involve division by sin or cos, so no need to
% special cases to avoid division by zero

% exp(i*phi),exp(-i*phi) are useful for both partial derivatives
expplus = exp(i*phi);
expminus = exp(-i*phi);

% theta derivative
% d/dtheta Y(n,m) = 1/2 exp(-i phi) sqrt((n-m)(n+m+1)) Y(n,m+1)
%                 - 1/2 exp(i phi) sqrt((n-m+1)(n+m)) Y(n,m-1)

ymplus = spharm(n,m+1,theta,phi);
ymminus = spharm(n,m-1,theta,phi);

Ytheta = sqrt((n-m+1)*(n+m))/2 * expplus .* ymminus ...
         - sqrt((n-m)*(n+m+1))/2 * expminus .* ymplus;

% phi derivative - actually 1/sin(theta) * d/dphi Y(n,m)
% Note that this is just i*m/sin(theta) * Y(n,m), but we use a
% recursion relation to avoid divide-by-zero trauma.
% 1/sin(theta) d/dphi Y(n,m) = 
% i/2 * [ exp(-i phi) sqrt((2n+1)(n+m+1)(n+m+2)/(2n+3)) Y(n+1,m+1)
%     + exp(i phi) sqrt((2n+1)(n-m+1)(n-m+2)/(2n+3)) Y(n+1,m-1) ]

ymplus = spharm(n+1,m+1,theta,phi);
ymminus = spharm(n+1,m-1,theta,phi);

Yphi = i/2 * sqrt((2*n+1)/(2*n+3)) * ...
   ( sqrt((n+m+1)*(n+m+2)) * expminus .* ymplus ...
   + sqrt((n-m+1)*(n-m+2)) * expplus .* ymminus );

return
