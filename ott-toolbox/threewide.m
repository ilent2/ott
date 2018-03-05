function wide_vector = threewide( a )% threewide.m - converts an input vector (either row or column%               vector) into a column vector repeated in%               three columns.% Usage:% wide_vector = threewide(original_vector);%% You might find this useful for multiplying a vector of scalars% with a column vector of 3-vectors.%% This file is part of the package Optical tweezers toolbox 1.0.1
% Copyright 2006-2007 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.htmla = a(:);wide_vector = [ a a a ];return