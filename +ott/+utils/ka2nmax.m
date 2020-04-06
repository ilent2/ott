function Nmax = ka2nmax(ka)
% Finds a reasonable Nmax to truncate at for given size parameter
%
% Usage
%   Nmax = ka2nmax(ka) calculates reasonable maximum order, Nmax, to
%   truncate beam beam coefficients/T-matrix at for a given size parameter.
%
% Returns :math:`Nmax = |ka| + 3 (|ka|)^(1/3)`

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ka = abs(ka);
Nmax = ka + 3 * ka.^(1/3);
Nmax = ceil(Nmax);
