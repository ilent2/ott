function x = isodd(input)
% ISODD determines if an integer is odd
% Outputs a matrix the same size as input with
% 1 for odd and 0 for even entries.
%
% Warning: Plays up if the the integer is of the order 10^16

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

x = logical(mod(input,2));
