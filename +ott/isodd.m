function x = isodd(value)
% Determines if an integer is odd
% Input: isodd(integer)
% Output: '1' if odd and '0' if even
%
% Warning: Plays up if the the integer is of the order 10^16
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

x = mod(value,2);

