function Lpl = laguerre(p,l,x)
% laguerre.m - associated Laguerre function
%
% Usage:
% L = laguerre(p,l,x)
%
% p and l must be integer scalars greater than zero
%
% Warning: this is a naive direct calculation, so might be slow or unstable
% for large p and/or l.
%
% This file is part of the package Optical tweezers toolbox 1.2
% Copyright 2006-2012 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

Lpl = nchoosek(p+l,p) * ones(size(x)); 

for m = 1:p
    Lpl = Lpl + (-1)^m/factorial(m) * nchoosek(p+l,p-m) * x.^m;
end

return
