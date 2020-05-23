classdef Gaussian < ott.beam.properties.Beam
% Properties of a paraxial Gaussian beam.
% Inherits from :class:`ott.beam.properties.Beam`.
%
% Properties
%   - waist         -- Beam waist radius
%   - power         -- Beam power

% Copyright 2020 Isaac Lenton
% This file is part of OTT, see LICENSE.md for information about
% using/distributing this file.

  properties
    waist          % Beam waist radius
  end

  methods % Getters/setters
    function beam = set.waist(beam, val)
      assert(isnumeric(val) && isscalar(val), ...
        'waist must be numeric scalar');
      beam.waist = val;
    end
  end
end
