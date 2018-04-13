function [hn,dhn] = sbesselh1(n,kr)
% SBESSELH1 spherical hankel function hn(kr) of the first kind,
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr).
%
% hn = SBESSELH1(n,z) calculates spherical hankel function of first kind.
%
% [hn,dzhn] = SBESSELH1(n,z) additionally, calculates the derivatives
% of the appropriate Ricatti-Bessel function divided by z.
%
% See also besselj and bessely.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:sbesselh1:move', ...
    'This file will move to ott.utils.sbesselh1');
ott_warning('internal');

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n-1];
end

[n,kr]=meshgrid(n,kr);

[hn] = besselh(n+1/2,1,kr);

hn = sqrt(pi./(2*kr)) .* (hn);

if nargout==2
    dhn=hn(1:end,end/2+1:end)-n(1:end,1:end/2)./kr(1:end,1:end/2) .* hn(1:end,1:end/2);
    hn=hn(1:end,1:end/2);
end

ott_warning('external');

return
