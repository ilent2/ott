function x = iseven(value)
% Determines if an integer is even
% Input: iseven(integer)
% Output: '1' if even and '0' if odd
%
% Warning: Plays up if the the integer is of the order 10^16
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

x = mod(mod(value,2)+1,2);

