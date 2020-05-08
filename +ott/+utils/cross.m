function M = cross(a, b)
% Calculate cross produce of two matrices
%
% Usage
%   M = cross(a, b)
%
% If either of a or b is a vector, replicates it first.

% TODO: There is another functionality we might want to implement too

assert(size(a, 1) == 3, 'num rows of a must be 3');
assert(size(b, 1) == 3, 'num rows of b must be 3');

assert((ndims(a) == ndims(b) && all(size(a) == size(b))) ...
    || isvector(a) || isvector(b), ...
    'a and b must match size or one must be a vector');

if isvector(a) && ~isvector(b)
  sz = size(b);
  a = repmat(a, 1, sz(2:end));
elseif isvector(b) && ~isvector(a)
  sz = size(a);
  b = repmat(b, [1, sz(2:end)]);
end

M = cross(a, b);

