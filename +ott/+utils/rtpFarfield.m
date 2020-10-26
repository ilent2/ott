function [tp, rtp] = rtpFarfield(rtp)
% Convert between 3xN and 2xN spherical coordinates with/without radius
%
% Usage
%   [tp, rtp] = rtpFarfield(input)
%
% Parameters
%   - input -- (2xN or 3xN numeric) Far-field theta/phi coordinates
%     with format [radius; theta; phi] or [theta; phi].
%
% Outputs
%   - tp -- (2xN numeric) 2xN formatted version with theta/phi.
%   - rtp -- (3xN numeric) 3xN formatted version (radius = 1).

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

assert(isnumeric(rtp) && ismatrix(rtp) ...
    && any(size(rtp, 1) == [2, 3]), ...
    'input must be 2xN or 3xN numeric matrix');

if size(rtp, 1) == 2
  tp = rtp;
else
  tp = rtp(2:3, :);
end

rtp = [ones(1, size(rtp, 2)); tp];

end

