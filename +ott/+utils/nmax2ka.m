function ka = nmax2ka(Nmax)
% Finds size parameter ka corresponding to Nmax.
%
% Truncation order is given by::
%
%   Nmax = ka + 3 (ka)^(1/3)
%
% Usage
%   ka = NMAX2KA(Nmax) finds size parameter for maximum order, Nmax,
%   which spherical expansions are truncated.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

assert(isnumeric(Nmax) && all(Nmax(:) >= 0) ...
  && all(round(Nmax(:)) == Nmax(:)), ...
  'Nmax must be positive numeric integers');

ka = zeros(size(Nmax));

for ii=1:length(Nmax)
    kas = roots([1 (-3*Nmax(ii)) (27+3*Nmax(ii).^2) (-Nmax(ii).^3)]);
    ka(ii) = kas(3);
end

ka(Nmax == 0) = 0;
