function [Aperp,Apar] = perpcomponent( A, n )
% PERPCOMPONENT finds perpendicular (and optionally) parallel
% components of a vector relative to a reference vector.
%
% perp_component = PERPCOMPONENT(A,n) calculates perpendicular component
% of row vector A or Nx3 matrix A of vectors.  n is the reference vector.
%
% [perp_component,parallel_component] = PERPCOMPONENT(A,n) calculates
% perpendicular and parallel components.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

import ott.utils.*

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
