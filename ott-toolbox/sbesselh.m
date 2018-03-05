function [hn,ierr] = sbesselh(n,htype,kr)
% sbesselh - spherical hankel function hn(kr)
%
% Usage:
% hn = sbesselh(n,htype,kr)
%
% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
%
% See besselh for more details
%
% This file is part of the package Optical tweezers toolbox 1.2
% Copyright 2006-2012 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html

kr=kr(:);
n=n(:);
[hn,ierr] = besselh(n'+1/2,htype,kr);
kr=repmat(kr,[1 length(n)]);
hn = sqrt(pi./(2*kr)) .* hn;

return
