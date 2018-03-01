function [T,T2,a,b] = tmatrix_mie(Nmax,k_medium,k_particle,radius)
% tmatrix_mie.m : Mie scattering and internal coefficients for uniform
%                 sphere, arranged as a sparse T-matrix.
%
% Usage:
% [T,T2] = tmatrix_mie(Nmax,k_medium,k_particle,radius)
% T is the T-matrix of scattered modes. T2 is the T-matrix of internal
% modes
%
% PACKAGE INFO

n=[1:Nmax];

m = k_particle/k_medium;

r0 = k_medium * radius;
r1 = k_particle * radius;

indexing=combined_index(1:Nmax^2+2*Nmax)';

j0 = (sbesselj(n,r0)).';
j1 = (sbesselj(n,r1)).';
h0 = (sbesselh1(n,r0)).';
j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

b = -( j1d.*j0 - m*j0d.*j1 ) ./ ( j1d.*h0 - m*h0d.*j1 );
a = -( j0d.*j1 - m*j1d.*j0 ) ./ ( h0d.*j1 - m*j1d.*h0 );

T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)],[a(indexing);b(indexing)]);

if nargout>1
    d = -( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
    c = -( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
    
    T2=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)],[c(indexing);d(indexing)]);
    
end