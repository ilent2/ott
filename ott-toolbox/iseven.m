function x = iseven(value)
% Determines if an integer is even
% Input: iseven(integer)
% Output: '1' if even and '0' if odd
%
% Warning: Plays up if the the integer is of the order 10^16
%
% PACKAGE_INFO

x = mod(mod(value,2)+1,2);

