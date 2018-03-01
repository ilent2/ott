function [hn] = sbesselh2(n,kr)
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

kr=kr(:);
n=n(:);
[hn] = besselh(n'+1/2,2,kr);
kr=repmat(kr,[1 length(n)]);
hn = sqrt(pi./(2*kr)) .* hn;

return
