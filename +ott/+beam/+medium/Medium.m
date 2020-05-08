classdef (Abstract) Medium
% Abstract base class for describing optical media
%
% Abstract properties
%   - permittivity  -- Permittivity of medium
%   - permeability  -- Permeability of medium
%   - speed         -- Wave speed in medium
%   - index         -- Refractive index of medium
%
% See also :class:`ott.beam.medium.Vacuum`, :class:`ott.beam.medium.Material`

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file

  properties (Abstract)
    permittivity
    permeability
    speed
    index
  end
end
