function Nmax = ka2nmax(ka)
% ka2nmax.m - Finds a reasonable maximum order to truncate at given
%             a size parameter ka
%
% Returns Nmax = ka + 3 (ka)^(1/3)
%
% PACKAGE INFO

Nmax = ka + 3 * ka.^(1/3);
Nmax = ceil(Nmax);

return
