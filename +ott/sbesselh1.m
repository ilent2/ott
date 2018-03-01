function [hn] = sbesselh1(n,kr)
% sbesselh1 - spherical hankel function hn(kr) of the first kind
%
% Usage:
% hn = sbesselh1(n,kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% See besselh for more details
%
% PACKAGE INFO

kr=kr(:);
n=n(:);
[hn] = besselh(n'+1/2,1,kr);
kr=repmat(kr,[1 length(n)]);
hn = sqrt(pi./(2*kr)) .* hn;

return
