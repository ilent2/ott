function rtp = sanitiseRtp(rtp)
% Convert spherical coordinates into range used by toolbox.
%
% The toolbox uses the convention ``r >= 0; 0 <= t <= pi; 0 <= p <= 2pi``.
%
% Usage
%   rtp = ott.utils.sanitiseRtp(rtp)

% Copyright 2020 Isaac Lenton
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% There are probably more optimal ways to do this conversion, but
% the following works fine
rtp = ott.utils.xyz2rtp(ott.utils.rtp2xyz(rtp));

