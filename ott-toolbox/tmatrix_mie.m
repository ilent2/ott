function T = tmatrix_mie(Nmax,k_medium,k_particle,radius)
% tmatrix_mie.m
% Mie scattering coefficients for uniform sphere, arranged as a
% sparse T-matrix
%
% Usage:
% T = tmatrix_mie(Nmax,k_medium,k_particle,radius)
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

m = k_particle/k_medium;
Z = 1/m;

r0 = k_medium * radius;
r1 = k_particle * radius;

t1 = [];
t2 = [];

for n = 1:Nmax

    j0 = sbesselj(n,r0);
    j1 = sbesselj(n,r1);
    j0d = sbesselj(n-1,r0) - n*sbesselj(n,r0)/r0;
    j1d = sbesselj(n-1,r1) - n*sbesselj(n,r1)/r1;
    h0 = sbesselh1(n,r0);
    h0d = sbesselh1(n-1,r0) - n*sbesselh1(n,r0)/r0;
    a = ( Z*j0d*j1 - j1d*j0 ) / ( j1d*h0 - Z*h0d*j1 );
    b = ( j0d*j1 - Z*j1d*j0 ) / ( Z*j1d*h0 - h0d*j1 );

    t1 = [ t1; a*ones(2*n+1,1) ];
    t2 = [ t2; b*ones(2*n+1,1) ];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spdiags didn't work for Adrian - looks like Matlab version
% incompatibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tsize = combined_index(Nmax,Nmax) * 2;

%T = spdiags([t1; t2],0,sparse([],[],[],tsize,tsize,tsize));

T = sparse(1:tsize,1:tsize,[t1; t2]);

return
