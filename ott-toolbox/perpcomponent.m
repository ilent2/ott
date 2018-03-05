function [Aperp,Apar] = perpcomponent( A, n )
% perpcomponent.m - finds perpendicular (and optionally) parallel
%       components of a vector relative to a reference vector.
%
% Usage:
% perp_component = perpcomponent(A,n)
% [perp_component,parallel_component] = perpcomponent(A,n)
%
% where n is the reference vector.
% The columns of A and n must be the vector components; a number
% of rows can be used to process many vectors at once.
%
% This file is part of the package Optical tweezers toolbox 1.2
% Copyright 2006-2012 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

% Make n a unit vector
nmag = dot(n,n,2);
nmag3 = threewide(nmag); % Only works for 3 component vectors
n = n./nmag3;

% The order in the dot product matters: the complex conjugate of
% the first term is used. Since we might want complex A ...
Aparmag = dot(n,A,2);
Aparmag3 = threewide(Aparmag);
Apar = Aparmag3 .* n;
Aperp = A - Apar;

return

