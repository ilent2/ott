function rcross = crossProductMatrix(rvec)
% Calculate the cross product matrix
%
% Calculates a 3x3 matrix for applying the cross product::
%
%   a = r \times p
%
% for a fixed vector r.
%
% Usage
%   rcross = crossProductMatrix(rvec)
%
% Returns
%   - rcross (3x3xN numeric) -- Array of 3x3 cross product matrices that can
%     be applied to a 3x1 vector.
%
% Parameters
%   - rvec (3xN numeric) -- The r component of the cross product.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

assert(ismatrix(rvec) && size(rvec, 1) == 3, ...
    'rvec must be 3xN matrix');

rcross = zeros(3, 3, size(rvec, 2));
rcross(1, 2, :) = -rvec(3, :);
rcross(1, 3, :) = rvec(2, :);
rcross(2, 3, :) = -rvec(1, :);
rcross(2, 1, :) = rvec(3, :);
rcross(3, 1, :) = -rvec(2, :);
rcross(3, 2, :) = rvec(1, :);

