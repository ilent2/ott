function [hn,dhn] = sbesselh(n,htype,kr)
% SBESSELH spherical hankel function hn(kr) of the first or 
% second kind hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% hn = SBESSELH(n,htype,z) computes the spherical hankel function
% of degree, n, of kind, htype, and argument, z. 
%
% [hn,dzhn] = SBESSELH(n,htype,z) additionally, calculates the 
% derivative of the appropriate Ricatti-Bessel function divided 
% by z.
%
% See also besselj and bessely.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott.warning('internal');

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n-1];
end

[n,kr]=meshgrid(n,kr);

[hn] = besselh(n+1/2,htype,kr);

hn = sqrt(pi./(2*kr)) .* (hn);

if nargout==2
    dhn=hn(1:end,end/2+1:end)-n(1:end,1:end/2)./kr(1:end,1:end/2) .* hn(1:end,1:end/2);
    hn=hn(1:end,1:end/2);
end

ott.warning('external');
