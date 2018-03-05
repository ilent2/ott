function [hn,ierr] = sbesselh2(n,kr)% sbesselh2 - spherical hankel function hn(kr) of the second kind%% Usage:% hn = sbesselh2(n,kr)%% hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)%% See besselh for more details%% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html[hn,ierr] = besselh(n+1/2,2,kr);hn = sqrt(pi./(2*kr)) .* hn;return
