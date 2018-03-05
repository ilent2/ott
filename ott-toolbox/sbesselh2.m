function [hn,ierr] = sbesselh2(n,kr)
% sbesselh2 - spherical hankel function hn(kr) of the second kind
%
% Usage:
% hn = sbesselh2(n,kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% See besselh for more details
%
% PACKAGE INFO

[hn,ierr] = besselh(n+1/2,2,kr);
hn = sqrt(pi./(2*kr)) .* hn;

return

