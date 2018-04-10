function [theta,phi] = angulargrid(ntheta,nphi,behaviour)
%ANGULARGRID makes a angular grid of points over a sphere
%
% [theta,phi] = ANGULARGRID(N) generates two column N^2-by-1 matrices
% with theta (polar) and phi (azimuthal) angle pairs for N discrete
% evenly spaced polar and azimuthal angles.
%
% [theta,phi] = ANGULARGRID(ntheta, nphi) specifies the number of
% evenly spaced points to use in the theta and phi direction.
%
% [theta,phi] = ANGULARGRID(..., behaviour) uses behaviour to control
% the output type:
%
% behaviour | output
% -----------------------------------------------
%     0     | column vectors of all points
%     1     | vectors of all theta and phi values 
%     2     | ntheta x nphi matrix of all points
%
% Note that the output data values are the same for
% behaviours 0 and 2; they're just arranged differently.
% To convert from one format to another:
% 2 -> 0: theta = theta(:); phi = phi(:);
% 0 -> 2: theta = reshape(theta,ntheta,nphi);
%         phi = reshape(phi,ntheta,nphi);
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

ott_warning('ott:angulargrid:move', ...
    'This file will move to ott.utils.angulargrid');

switch nargin
case 1
   nphi = ntheta;
   behaviour = 0;
case 2
   behaviour = 0;
otherwise
   % Everything has been specified
end

% theta goes from 0 to pi - we avoid the endpoints
% since they are mathematically troublesome
theta = ((1:ntheta)-0.5)/ntheta * pi;
% phi goes from 0 to 2*pi, 0 = 2*pi
phi = ((1:nphi)-1)/nphi * 2*pi;

if behaviour == 1
   % All we need are the vectors of the theta and phi values
   return
end

theta = theta.' * ones(1,nphi);
phi = ones(ntheta,1) * phi;

if behaviour == 2
   return
end

theta = theta(:);
phi = phi(:);

return

