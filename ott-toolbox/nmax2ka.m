function ka = nmax2ka(Nmax)
% nmax2ka.m - Finds size parameter ka corresponding to maximum order at
%             which spherical expansions are truncated
%
% Truncation order is given by Nmax = ka + 3 (ka)^(1/3)
%
% This file is part of the package Optical tweezers toolbox 1.2
% Copyright 2006-2012 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

for ii=1:length(Nmax)
    kas = roots([Nmax(ii) (-3*Nmax(ii)) (27+3*Nmax(ii).^2) (-Nmax(ii).^3)]);
    ka(ii) = kas(3);
end

return
