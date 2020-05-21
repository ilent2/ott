function [hn,dhn] = sbesselh2(n,kr, varargin)
% SBESSELH2 spherical hankel function hn(kr) of the second kind
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% hn = SBESSELH2(n,z) calculates the hankel function of the second kind.
%
% [hn,dzhn] = SBESSELH2(n,z) additionally, calculates the derivative
% of the appropriate Ricatti-Bessel function divided by z.
%
% See also besselj and bessely.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott.warning('internal');

p = inputParser;
p.addParameter('method', []);
p.parse(varargin{:});

method = p.Results.method;
if isempty(method)
  if any(abs(kr) < max(abs(n)))
    method = 'besselj+bessely';
  else
    method = 'besselh';
  end
end

kr=kr(:);
n=n(:);

if nargout==2
    n=[n;n-1];
end

[n,kr]=meshgrid(n,kr);

switch method

  case 'besselj+bessely'
    bj = besselj(n+1/2,kr);
    by = bessely(n+1/2,kr);
    hn = bj - 1i*by;

  case 'besselh'
    % This doesn't work too well for kr < n
    [hn] = besselh(n+1/2,2,kr);

  otherwise
    error('Unknown method');
end

hn = sqrt(pi./(2*kr)) .* (hn);

if nargout==2
    dhn=hn(1:end,end/2+1:end)-n(1:end,1:end/2)./kr(1:end,1:end/2) .* hn(1:end,1:end/2);
    hn=hn(1:end,1:end/2);
end

ott.warning('external');
