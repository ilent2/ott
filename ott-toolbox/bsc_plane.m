function [nn,mm,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
% Finds VSWF representation of plane wave
% [n,m,a,b] = bsc_plane(nmax,lambda,theta,phi,Etheta,Ephi)
%
% theta and phi (in radians) give the direction of propagation
% of the plane wave. +z direction is theta = 0, phi = any value
%
% PACKAGE_INFO

E = [ 0 Etheta Ephi ];

[nn,mm]=combined_index([1:nmax*(nmax+2)].');

a = zeros(size(nn));
b = zeros(size(nn));

ii=0;
for n = 1:nmax
    
    iter=2*n+1;
    
    Nn = sqrt(1/(n*(n+1)));
    
    [B,C,P] = vsh(n,[-n:n],theta,phi);
    
    Creorg=reshape(C,[iter,3]);
    Breorg=reshape(B,[iter,3]);
    
    Et=repmat(E,[iter,1]);
    
    % dot product function takes complex conjugate of first term
    a(ii+1:ii+iter) = 4*pi * 1i^n * Nn * dot(Creorg,Et,2);
    b(ii+1:ii+iter) = 4*pi * 1i^(n-1) * Nn * dot(Breorg,Et,2);
    
    ii=ii+iter;
    
end

return
