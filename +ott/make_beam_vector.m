function [newa,newb,newn,newm] = make_beam_vector(olda,oldb,n,m,Nmax)
% make_beam_vector.m
%
% Convert the n,m,a,b as output by the bsc_* functions to sparse vector
% a and b with standard packing
%
% Usage:
% [a2,b2] = make_beam_vector(a1,b1,n,m);
% [a2,b2] = make_beam_vector(a1,b1,n,m,Nmax);
%
% Optional (will not calculate if n2, m2 not requested):
% [a2,b2,n2,m2] = make_beam_vector(a1,b1,n,m);
% [a2,b2,n2,m2] = make_beam_vector(a1,b1,n,m,Nmax);
%
% PACKAGE INFO

if nargin < 5
    Nmax = max(n);
end

total_orders = combined_index(Nmax,Nmax);
ci = combined_index(n,m);

newa = sparse(ci,1,olda,total_orders,1);
newb = sparse(ci,1,oldb,total_orders,1);

if nargout>2
    [newn,newm]=combined_index(1:Nmax^2+2*Nmax);
    
    newn=newn.';
    newm=newm.';
end

return
