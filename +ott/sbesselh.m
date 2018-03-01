function [hn] = sbesselh(n,htype,kr)
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

kr=kr(:);
n=n(:);
[hn] = besselh(n'+1/2,htype,kr);
kr=repmat(kr,[1 length(n)]);
hn = sqrt(pi./(2*kr)) .* hn;

return
