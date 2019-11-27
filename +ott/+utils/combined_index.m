function [out1,out2] = combined_index(in1,in2)
%COMBINED_INDEX translates between (n,m) and combined index
% Mode indices and combined index are related by: ci = n * (n+1) + m.
%
% [n,m] = COMBINED_INDEX(ci) calculates (n,m) from the combined index.
%
% ci = COMBINED_INDEX(n,m) calculates the combined index from mode indices.
%
% length = COMBINED_INDEX(Nmax, Nmax) calculates length of the beam vectors.
%
% Nmax = COMBINED_INDEX(length) calculates Nmax from length of beam vectors.

% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Sanity check
if nargin == 1
   out1 = floor(sqrt(in1));
   out2 = in1 - out1.^2 - out1;
elseif nargin == 2
   out1 = in1 .* (in1 + 1) + in2;
else
   error('Bad number of input/output arguments');
end
