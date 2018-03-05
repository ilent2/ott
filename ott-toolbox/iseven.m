function x = iseven(input)
% ISEVEN determines if an integer is even
% Outputs a matrix of the same size as input with 1
% for even and 0 for odd entries.
%
% Warning: Plays up if the the integer is of the order 10^16
%
% PACKAGE_INFO

warning('This function will be moved to ott.utils.iseven');

x = mod(mod(input,2)+1,2);

