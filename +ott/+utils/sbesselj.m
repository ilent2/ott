function [jn,djn] = sbesselj(n,kr)
% SBESSELJ spherical bessel function jn(kr)
% jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
%
% jn = SBESSEL(n,z) calculates the spherical bessel function.
%
% [jn,dzjn] = sbessel(n,z) additionally, calculates the derivative
% of the appropriate Ricatti-Bessel function divided by z.
%
% When kr is less than 1e-15, uses a small value approximation invovling
% n double factorial.
%
% See also besselj.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.*

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
    jn = kr.^n ./ prod(1:2:(2*n+1));
elseif length(kr) == 1 && ~(abs(kr) < 1e-15)
    jn = sqrt(pi./(2*kr)) .* jn;
elseif length(n) == 1
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n ./ prod(1:2:(2*n+1));
else % both n and kr are vectors
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
      
    % Calculate n!! (previous definition used max(n) only)
    npp = 1:2:(2*max(n(small_args))+1);
    npp = repmat(npp, numel(n(small_args)), 1);
    npp(npp > reshape(2*n(small_args)+1, [], 1)) = 1;
    npp = reshape(prod(npp, 2), size(n(small_args)));
      
    jn(small_args) = kr(small_args).^n(small_args) ./ npp;
end

if nargout==2
    djn=jn(:,end/2+1:end)-n(:,1:end/2)./kr(:,1:end/2) .* jn(:,1:end/2);
    jn=jn(:,1:end/2);
end

