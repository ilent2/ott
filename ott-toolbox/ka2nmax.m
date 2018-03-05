function Nmax = ka2nmax(ka)
% ka2nmax.m - Finds a reasonable maximum order to truncate at given
%             a size parameter ka
%
% Returns Nmax = ka + 3 (ka)^(1/3)
%
% This file is part of the package Optical tweezers toolbox 1.0.1
% Copyright 2006-2007 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

Nmax = ka + 3 * ka^(1/3);
Nmax = ceil(Nmax);

return
