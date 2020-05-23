function M = cart2sph_mat(theta, phi)
% Generate Cartesian to spherical coordinate conversion matrix
%
% Usage
%   M = cart2sph_mat(theta, phi)
%
% See also :func:`ott.utils.sph2cart_mat`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

M = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
     cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
     -sin(phi) cos(phi) 0];

end
