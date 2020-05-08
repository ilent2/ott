function M = dot(a, b)
% Calculate dot product of two matrices
%
% Usage
%   M = dot(a, b)
%
% If either of a or b is a vector, replicates it first.

% TODO: There is another functionality we might want to implement too

assert(size(a, 1) == 3, 'num rows of a must be 3');
assert(size(b, 1) == 3, 'num rows of b must be 3');

assert((ndims(a) == ndims(b) && all(size(a) == size(b))) ...
    || isvector(a) || isvector(b), ...
    'a and b must match size or one must be a vector');

M = sum(a .* b, 1);

