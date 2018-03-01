function Lpl = laguerre(p,l,x)
% laguerre.m - associated Laguerre function
%
% Usage:
% L = laguerre(p,l,x)
%
% p and l must be integer scalars greater than zero
%
% Warning: this is a naive direct calculation, so might be slow or unstable
% for large p and/or l.
%
% PACKAGE INFO

Lpl = nchoosek(p+l,p) * ones(size(x)); 

for m = 1:p
    Lpl = Lpl + (-1)^m/factorial(m) * nchoosek(p+l,p-m) * x.^m;
end

return
