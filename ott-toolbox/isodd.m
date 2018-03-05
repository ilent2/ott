function x = isodd(value)
% Determines if an integer is odd
% Input: isodd(integer)
% Output: '1' if odd and '0' if even
%
% Warning: Plays up if the the integer is of the order 10^16
%
% PACKAGE_INFO

x = mod(value,2);

