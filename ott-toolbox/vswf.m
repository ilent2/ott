function [M,N,M2,N2,M3,N3] = vswf(n,m,kr,theta,phi,htype)
% vswf.m : vector spherical harmonics
%
% Usage:
% [M,N] = vswf(n,m,kr,theta,phi,type)
% or
% [M1,N1,M2,N2,M3,N3] = vswf(n,m,kr,theta,phi)
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
% The three components of each vector are [r,theta,phi]
%
% "Out of range" n and m result in return of [0 0 0]
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

% Check input vectors
% These must all be of equal length if non-scalar
% and for good measure, we expand any scalar ones
% to match the others in length

% Convert all to column vectors
kr = kr(:);
theta = theta(:);
phi = phi(:);

% Check the lengths
[kr,theta,phi] = matchsize(kr,theta,phi);

if nargin < 6
   htype = 0;
end

[B,C,P] = vsh(n,m,theta,phi);
if n > 0
   Nn = sqrt(1/(n*(n+1)));
else
   Nn = 0;
end

switch(htype)
case 1,
   hn = threewide(sbesselh1(n,kr));
   hn1 = threewide(sbesselh1(n-1,kr));
   kr3 = threewide(kr);
   M = Nn * hn .* C;
   N = Nn * ( n*(n+1)./kr3 .* hn .* P + ( hn1 - n./kr3 .* hn ) .* B );
   M2 = 0; N2 = 0; M3 = 0; N3 = 0;
case 2,
   hn = threewide(sbesselh2(n,kr));
   hn1 = threewide(sbesselh2(n-1,kr));
   kr3 = threewide(kr);
   M = Nn * hn .* C;
   N = Nn * ( n*(n+1)./kr3 .* hn .* P + ( hn1 - n./kr3 .* hn ) .* B );
   M2 = 0; N2 = 0; M3 = 0; N3 = 0;
case 3,
   jn = threewide(sbesselj(n,kr));
   jn1 = threewide(sbesselj(n-1,kr));
   kr3 = threewide(kr);
   M = Nn * jn .* C;
   N = Nn * ( n*(n+1)./kr3 .* jn .* P + ( jn1 - n./kr3 .* jn ) .* B );
   M2 = 0; N2 = 0; M3 = 0; N3 = 0;
otherwise,
   hn1 = threewide(sbesselh1(n,kr));
   hn2 = threewide(sbesselh2(n,kr));
   jn = threewide(sbesselj(n,kr));
   hn11 = threewide(sbesselh1(n-1,kr));
   hn21 = threewide(sbesselh2(n-1,kr));
   jn1 = threewide(sbesselj(n-1,kr));
   kr3 = threewide(kr);
   M = Nn * hn1 .* C;
   N = Nn * ( n*(n+1)./kr3 .* hn1 .* P + ( hn11 - n./kr3 .* hn1 ) .* B );
   M2 = Nn * hn2 .* C;
   N2 = Nn * ( n*(n+1)./kr3 .* hn2 .* P + ( hn21 - n./kr3 .* hn2 ) .* B );
   M3 = Nn * jn .* C;
   N3 = Nn * ( n*(n+1)./kr3 .* jn .* P + ( jn1 - n./kr3 .* jn ) .* B );
end

return

