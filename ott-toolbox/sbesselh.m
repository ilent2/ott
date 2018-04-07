function [hn,dhn] = sbesselh(n,htype,kr)
% sbesselh - spherical hankel function hn^(htype)(kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% Usage:
%
% hn = sbesselh(n,z);
% OR
% [hn,dzhn] = sbesselh(n,z);
%
% where dzhn is the derivative of the appropriate Ricatti-Bessel function
% divided by z.
%
% See besselj and bessely for more details
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n-1];
end

[n,kr]=meshgrid(n,kr);

[jn] = besselj(n+1/2,kr);
[yn] = bessely(n+1/2,kr);

hn = sqrt(pi./(2*kr)) .* (jn+1i*(2*(htype==1)-1)*yn);

if nargout==2
    dhn=hn(1:end,end/2+1:end)-n(1:end,1:end/2)./kr(1:end,1:end/2) .* hn(1:end,1:end/2);
    hn=hn(1:end,1:end/2);
end



return
