function Nmax = ka2nmax(ka)
% ka2nmax.m - Finds a reasonable maximum order to truncate at given
%             a size parameter ka
%
% Returns Nmax = ka + 3 (ka)^(1/3)
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

Nmax = ka + 3 * ka.^(1/3);
Nmax = ceil(Nmax);

return
