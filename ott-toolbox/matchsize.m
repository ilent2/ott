function [A,B,C] = matchsize(A,B,C)
% MATCHSIZE checks that all vector inputs have the same number of rows
%
% [A,B] = MATCHSIZE(A,B) checks inputs have same number of rows,
% and expands single-row inputs by repetition to match the input row number.
%
% [A,B,C] = matchsize(A,B,C) as above but for 3 inputs/outputs.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

An = size(A);
Bn = size(B);
if nargin >= 3
   Cn = size(C);
else
   Cn = [1 1];
   C = [0];
end

nmax = max(An(1),Bn(1));
nmax = max(nmax,Cn(1));

if An(1) < nmax
   if An(1) == 1
      A = ones(nmax,1) * A;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end
if Bn(1) < nmax
   if Bn(1) == 1
      B = ones(nmax,1) * B;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end
if Cn(1) < nmax & nargin >= 3
   if Cn(1) == 1
      C = ones(nmax,1) * C;
  else
      error(['Number of rows in inputs must be one or equal.']);
  end
end

return

