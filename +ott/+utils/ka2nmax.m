function Nmax = ka2nmax(ka)
% KA2NMAX finds a reasonable Nmax to truncate at for given size parameter
%
% Nmax = KA2NMAX(ka) calculates reasonable maximum order, Nmax, to
% truncate beam beam coefficients/T-matrix at for a given size parameter.
%
% Returns Nmax = ka + 3 (ka)^(1/3)
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

Nmax = ka + 3 * ka.^(1/3);
Nmax = ceil(Nmax);
