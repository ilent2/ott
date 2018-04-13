function [jn,djn] = sbesselj(n,kr)
% SBESSELJ spherical bessel function jn(kr)
% jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
%
% jn = SBESSEL(n,z) calculates the spherical bessel function.
%
% [jn,dzjn] = sbessel(n,z) additionally, calculates the derivative
% of the appropriate Ricatti-Bessel function divided by z.
%
% See also besselj.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:sbesselj:move', ...
    'This file will move to ott.utils.sbesselj');
ott_warning('internal');

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n-1];
end

[n,kr]=meshgrid(n,kr);
[jn] = besselj(n+1/2,kr);

small_args = find( abs(kr) < 1e-15 );
not_small_args = find( ~(abs(kr) < 1e-15) );

if length(kr) == 1 && abs(kr) < 1e-15
    jn = kr.^n ./ repmat(prod(1:2:(2*n+1)),[kr,1]);
elseif length(kr) == 1 && ~(abs(kr) < 1e-15)
    jn = sqrt(pi./(2*kr)) .* jn;
elseif length(n) == 1
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n ./ prod(1:2:(2*n+1));
else % both n and kr are vectors
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n(small_args) ./ ...
        prod(1:2:(2*n(small_args)+1));
end

if nargout==2
    djn=jn(1:end,end/2+1:end)-n(1:end,1:end/2)./kr(1:end,1:end/2) .* jn(1:end,1:end/2);
    jn=jn(1:end,1:end/2);
end

ott_warning('external');

return
