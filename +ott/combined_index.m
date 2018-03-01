function [out1,out2] = combined_index(in1,in2)
% combined_index.m - translates between (n,m) and combined index
%                    ci = n * (n+1) + m
%
% Usage:
% [n,m] = combined_index(ci);
% ci = combined_index(n,m);
%
% This file is part of the package Optical tweezers toolbox 1.3
% Copyright 2006-2013 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

% Sanity check
if nargin == 1
   out1 = floor(sqrt(in1));
   out2 = in1 - out1.^2 - out1;
elseif nargin == 2
   out1 = in1 .* (in1 + 1) + in2;
else
   error('Bad number of input/output arguments');
end

return
