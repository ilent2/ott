function [out1,out2] = combined_index(in1,in2)
% combined_index.m - translates between (n,m) and combined index
%                    ci = n * (n+1) + m
%
% Usage:
% [n,m] = combined_index(ci);
% ci = combined_index(n,m);
%
% PACKAGE_INFO

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
