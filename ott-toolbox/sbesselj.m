function [jn,ierr] = sbesselj(n,kr)
% sbesselj - spherical bessel function jn(kr)
%
% jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
%
% Usage: jn = sbessel(n,z); 
%
% See besselj for more details
%
% This file is part of the package Optical tweezers toolbox 1.2
% Copyright 2006-2012 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

kr=kr(:);
n=n(:);
[jn,ierr] = besselj(n'+1/2,kr);
[n,kr]=meshgrid(n,kr);

small_args = find( abs(kr) < 1e-15 );
not_small_args = find( ~(abs(kr) < 1e-15) );

if length(kr) == 1 & abs(kr) < 1e-15
    jn = kr.^n ./ prod(1:2:(2*n+1));
elseif length(kr) == 1 & ~(abs(kr) < 1e-15)
    jn = sqrt(pi./(2*kr)) .* jn;
elseif length(n) == 1
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n ./ prod(1:2:(2*n+1));
else % both n and kr are vectors
    jn(not_small_args) = ...
        sqrt(pi./(2*kr(not_small_args))) .* jn(not_small_args);
    jn(small_args) = kr(small_args).^n(small_args) ./ ...
        prod(1:2:(2*n(small_args)+1));
end

return
