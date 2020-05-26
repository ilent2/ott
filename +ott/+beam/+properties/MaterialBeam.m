classdef (Abstract) MaterialBeam < ott.beam.properties.Beam ...
    & ott.beam.properties.Material
% Base class for beams defining material properties.
% Inherits from :class:`Beam` and :class:`ott.beam.properties.MaterialBeam`.
%
% Properties
%   - omega       -- Optical frequency of beam.
%   - medium      -- Properties of optical medium.

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  methods (Static)
    function args = likeProperties(other, args)
      args = ott.beam.properties.Material.likeProperties(other, args);
      args = ott.beam.properties.Beam.likeProperties(other, args);
    end
  end
end
