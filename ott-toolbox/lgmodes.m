function A = lgmodes(r,phi,E,maxp,maxl)
% lgmodes -- Decomposition of paraxial beam into LG modes
%
% Usage:
% A = lgmodes(r,phi,E,maxp,maxl);
% where
% r is in units of the beam width,
% E is in arbitrary units,
% r,phi,E are either vectors of equal length, or r and phi are
% vectors and E is a matrix of size length(r) by length(phi)
% maxp, maxl are the maximum radial and azimuthal mode indices
% A is a matrix of the mode amplitudes with elements
% Apl (ie rows are radial modes, the first row is p = 0,
% columns are azimuthal modes, note that l = 0 is the central column.)
%
% WARNING: INCOMPLETE!
% r,phi,E MUST BE COLUMN VECTORS OF EQUAL LENGTH
%
% PACKAGE INFO

% Is E a vector or a matrix?
% For now, assume that r,phi,E are column vectors

number_of_modes = ( maxp + 1 ) * ( 2*maxl + 1 );
number_of_points = length(r);

coefficient_matrix = zeros(number_of_points,number_of_modes);

p = 0:maxp;
l = -maxl:maxl;

p = p(:) * ones(1,2*maxl+1);
l = ones(maxp+1,1) * l;

p = p(:);
l = l(:);

for n = 1:number_of_modes
    coefficient_matrix(:,n) = lgmode(p(n),l(n),r,phi);
end

A = coefficient_matrix\E;

%[p,l,A]

A = reshape(A,maxp+1,2*maxl+1);

return
