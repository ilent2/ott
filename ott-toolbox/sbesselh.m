function [hn,ierr] = sbesselh(n,htype,kr)
% sbesselh - spherical hankel function hn(kr)
%
% Usage:
% hn = sbesselh(n,htype,kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% See besselh for more details
%
% PACKAGE INFO

[hn,ierr] = besselh(n+1/2,htype,kr);
hn = sqrt(pi./(2*kr)) .* hn;

return

