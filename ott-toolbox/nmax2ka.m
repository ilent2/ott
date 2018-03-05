function ka = nmax2ka(Nmax)
% nmax2ka.m - Finds size parameter ka corresponding to maximum order at
%             which spherical expansions are truncated
%
% Truncation order is given by Nmax = ka + 3 (ka)^(1/3)
%
% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

kas = roots([1 (-3*Nmax) (27+3*Nmax^2) (-Nmax^3)]);
ka = kas(3);

return
