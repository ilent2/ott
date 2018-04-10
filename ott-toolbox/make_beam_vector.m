function [newa,newb,newn,newm] = make_beam_vector(olda,oldb,n,m,Nmax)
% MAKE_BEAM_VECTOR convert output of bsc_* functions to sparse vectors.
%
% [a2,b2,n2,m2] = make_beam_vector(a1,b1,n,m,Nmax) generates sparse
% beam vectors a2, b2 using standard packing.  The length of the beam
% vectors is set by Nmax.
%
% [a2,b2,n2,m2] = make_beam_vector(a1,b1,n,m) calculates Nmax from
% max(n).
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

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
