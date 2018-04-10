function wide_vector = threewide( a )
% THREEWIDE creates colum vector with input repeated in 3 columns
%   the function can take a column of row vector input, the output
%   will be a matrix with three columns.
%
% You might find this useful for multiplying a vector of scalars
% with a column vector of 3-vectors.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:threewide:move', ...
    'this function will move to ott.utils.threewide');

a = a(:);
wide_vector = [ a a a ];
